/**
 * @author Rushabh Tolia, Dan Rattanakornphan, Wesley Burnawan
 *
 * @section LICENSE
 * Copyright (c) 2012 Database Group, Computer Sciences Department, University of Wisconsin-Madison.
 */

#include "btree.h"
#include "filescan.h"
#include "exceptions/bad_index_info_exception.h"
#include "exceptions/bad_opcodes_exception.h"
#include "exceptions/bad_scanrange_exception.h"
#include "exceptions/no_such_key_found_exception.h"
#include "exceptions/scan_not_initialized_exception.h"
#include "exceptions/index_scan_completed_exception.h"
#include "exceptions/file_not_found_exception.h"
#include "exceptions/end_of_file_exception.h"


//#define DEBUG

namespace badgerdb
{
  PageId firstRootValue;
  const PageId INVALID_NUMBER = 0;

  // -----------------------------------------------------------------------------
  // BTreeIndex::BTreeIndex -- Constructor
  // -----------------------------------------------------------------------------

  /**
   * BTreeIndex Constructor. 
   * Check to see if the corresponding index file exists. If so, open the file.
   * If not, create it and insert entries for every tuple in the base relation using FileScan class.
   *
   * @param relationName        Name of file.
   * @param outIndexName        Return the name of index file.
   * @param bufMgrIn                        Buffer Manager Instance
   * @param attrByteOffset          Offset of attribute, over which index is to be built, in the record
   * @param attrType                        Datatype of attribute over which index is built
   */
BTreeIndex::BTreeIndex(const std::string & relationName,
        std::string & outIndexName,
        BufMgr *bufMgrIn,
        const int attrByteOffset,
        const Datatype attrType)
{
   nodeOccupancy = INTARRAYNONLEAFSIZE;
   leafOccupancy = INTARRAYLEAFSIZE;
   bufMgr = bufMgrIn;
   scanExecuting = false;

   std::ostringstream idxStr;
   idxStr << relationName << "." << attrByteOffset;
   outIndexName = idxStr.str();

   try {
	   Page *firstPageInBlob;
       file = new BlobFile(outIndexName, false);
       headerPageNum = file->getFirstPageNo();
       bufMgr->readPage(file, headerPageNum, firstPageInBlob);
       IndexMetaInfo *metaPage = (IndexMetaInfo *)firstPageInBlob;
       rootPageNum = metaPage->rootPageNo;

       if (relationName != metaPage->relationName || attrByteOffset != metaPage->attrByteOffset
           || attrType != metaPage->attrType) {
         throw BadIndexInfoException(outIndexName);
       }
       bufMgr->unPinPage(file, headerPageNum, false);
     } catch(FileNotFoundException e) {
       Page *rootPage;
       Page *firstPage;
	   file = new BlobFile(outIndexName, true);
       bufMgr->allocPage(file, rootPageNum, rootPage);
       bufMgr->allocPage(file, headerPageNum, firstPage);

       IndexMetaInfo *metaPage = (IndexMetaInfo *)firstPage;
       metaPage->rootPageNo = rootPageNum;
       metaPage->attrByteOffset = attrByteOffset;
       metaPage->attrType = attrType;
       relationName.copy(metaPage->relationName, 20, 0);
       
       LeafNodeInt *root = (LeafNodeInt *)rootPage;
       firstRootValue = rootPageNum;
	   root->rightSibPageNo = INVALID_NUMBER;

       bufMgr->unPinPage(file, rootPageNum, true);
	   bufMgr->unPinPage(file, headerPageNum, true);

	   RecordId recordID;
       FileScan scannedFile(relationName, bufMgr);

       try{
         for(;;) {
           scannedFile.scanNext(recordID);
           std::string currentRecord = scannedFile.getRecord();
           insertEntry(currentRecord.c_str() + attrByteOffset, recordID);
         }
       }
       catch(EndOfFileException e) {
         bufMgr->flushFile(file);
       }
     }
}

// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------

  /**
   * BTreeIndex Destructor. 
   * End any initialized scan, flush index file, after unpinning any pinned pages, from the buffer manager
   * and delete file instance thereby closing the index file.
   * Destructor should not throw any exceptions. All exceptions should be caught in here itself.
   **/
BTreeIndex::~BTreeIndex()
{
  scanExecuting = false;
  bufMgr->flushFile(BTreeIndex::file);
  delete file;
  file = nullptr;
}

/*
// The follwing method is used to help insert data into the btree by looking for the matching leaf node
// based on the internal nodes.
//
// @param nonLeafNode, the internal node being searched for
// @param insertedData, the data that is trying to be entered in.
*/
PageId BTreeIndex::findBelowNode(NonLeafNodeInt *nonLeafNode, int insertedData)
{
  int i = 0;
  if(insertedData < nonLeafNode->keyArray[i]) {
    return nonLeafNode->pageNoArray[i];
  }
  while (i < nodeOccupancy && nonLeafNode->pageNoArray[i+1] != INVALID_NUMBER && (nonLeafNode->keyArray[i] <= insertedData)) {
    i++;
  }
    return nonLeafNode->pageNoArray[i];
}

/*
// The follwing method is recursively called to insert data into the corresponding file to help insert nodes.
//
// @param key, the value of the key to find the data
// @param rid, the matching record id
*/
void BTreeIndex::insertData(Page *currentPage, PageId currentNodeID, bool isLeaf, const RIDKeyPair<int> insertedData, PageKeyPair<int> *&pushedUpEntry)
{
  PageId belowNodeID;
  Page *belowNode;	
  if (isLeaf) {
    LeafNodeInt *leafNode = (LeafNodeInt *)currentPage;
    if (leafNode->ridArray[leafOccupancy - 1].page_number == 0)
    {
      childEntry(leafNode, insertedData);
      bufMgr->unPinPage(file, currentNodeID, true);
      pushedUpEntry = nullptr;
    }
    else {
      childSplit(leafNode, currentNodeID, pushedUpEntry, insertedData);
    }
  }
  else{
    NonLeafNodeInt *nonLeafNode = (NonLeafNodeInt *)currentPage;
    belowNodeID = findBelowNode(nonLeafNode, insertedData.key);
    bufMgr->readPage(file, belowNodeID, belowNode);
    if (nonLeafNode->level == 1) {
      isLeaf = true;
    } else {
      isLeaf = false;
    }
    insertData(belowNode, belowNodeID, isLeaf, insertedData, pushedUpEntry);
    if (pushedUpEntry == nullptr) {
        bufMgr->unPinPage(file, currentNodeID, false);
        return;
    }
      if (nonLeafNode->pageNoArray[nodeOccupancy] == INVALID_NUMBER) {
        internalNodeEntry(nonLeafNode, pushedUpEntry);
        pushedUpEntry = nullptr;
        bufMgr->unPinPage(file, currentNodeID, true);
      }
      else
      {
        internalNodeSplit(nonLeafNode, currentNodeID, pushedUpEntry);
      }
  }
}

/*
// The follwing method is used to help split the internal nodes rather than the leaf nodes. 
//
// @param leftNode, the node that is getting ready to be split 
// @param leftNodeID, the record id of the node left to the current node
// @param pushedUpEntry, the node that get pushed up after following the split of an internal node
*/
void BTreeIndex::internalNodeSplit(NonLeafNodeInt *leftNode, PageId leftNodeID, PageKeyPair<int> *&pushedUpEntry)
{
  PageId rightNodeID;
  Page *rightPage;
  NonLeafNodeInt *rightNode;
  int tmpKey;
  bufMgr->allocPage(file, rightNodeID, rightPage);
  rightNode = (NonLeafNodeInt *)rightPage;
  int middleIndex = nodeOccupancy / 2;
  int pushupIndex = middleIndex;

  if (nodeOccupancy % 2 == 0 && leftNode->keyArray[middleIndex] > pushedUpEntry->key) {
    pushupIndex = middleIndex - 1;
  }  else {
   middleIndex++;
   }

  tmpKey = leftNode->keyArray[pushupIndex];
  for(int i = middleIndex; i < nodeOccupancy; i++) {
    rightNode->keyArray[i-middleIndex] = leftNode->keyArray[i];
    rightNode->pageNoArray[i-middleIndex] = leftNode->pageNoArray[i+1];
    leftNode->keyArray[i+1] = 0;
	leftNode->pageNoArray[i+1] = INVALID_NUMBER;
  }
  leftNode->keyArray[pushupIndex] = 0;
  leftNode->pageNoArray[pushupIndex] = INVALID_NUMBER;
  rightNode->level = leftNode->level;

  if (rightNode->keyArray[0] > pushedUpEntry->key) {
    internalNodeEntry(leftNode, pushedUpEntry);
  } else {
    internalNodeEntry(rightNode, pushedUpEntry);
  }
  pushedUpEntry->set(rightNodeID, tmpKey);
  bufMgr->unPinPage(file, rightNodeID, true);
  bufMgr->unPinPage(file, leftNodeID, true);

  if (rootPageNum == leftNodeID) {
    createNewRoot(leftNodeID, pushedUpEntry, false);
  }
}

/*
// The follwing method is used help create a new root following a splitting of the orirignal into a internal node. 
//
// @param formerRootID, the value of the original root id that was placed within the btree
// @param newRoot, the new root that will be added in place of the original one
// @param isLeaf, a boolean value for whether the node is a leaf or an internal node
*/
void BTreeIndex::createNewRoot(PageId formerRootID, PageKeyPair<int> *newRoot, bool isLeaf)
{
  PageId newRootNodeID;
  Page *newRootPage;
  Page *firstPageInBlob;
  bufMgr->allocPage(file, newRootNodeID, newRootPage);
  NonLeafNodeInt *newRootNode = (NonLeafNodeInt *)newRootPage;

  newRootNode->pageNoArray[0] = formerRootID;
  newRootNode->pageNoArray[1] = newRoot->pageNo;
  newRootNode->keyArray[0] = newRoot->key;

  if (isLeaf) {
    newRootNode->level = 1;
  } else {
    newRootNode->level = 0;
  }
  IndexMetaInfo *metaPage;
  bufMgr->readPage(file, headerPageNum, firstPageInBlob);
  metaPage = (IndexMetaInfo *)firstPageInBlob;
  metaPage->rootPageNo = newRootNodeID;
  rootPageNum = newRootNodeID;
  bufMgr->unPinPage(file, newRootNodeID, true);
  bufMgr->unPinPage(file, headerPageNum, true);
}

/*
// The following function is a helper to split a child node once full.
//
// @param leafNode, the present value of the leaf where the splitting must occur
// @param leafNodeID, the matching id of the leaf node needing to be split
// @param pushedUpEntry, the entry that will need to be pushed up accordingly
// @param dataEntry, the data value that is trying to get entered into the Btree
*/
void BTreeIndex::childSplit(LeafNodeInt *leafNode, PageId leafNodeID, PageKeyPair<int> *&pushedUpEntry, const RIDKeyPair<int> dataEntry)
{
  int middleIndex = leafOccupancy/2;
  PageId newLeafID;
  Page *newLeaf;
  LeafNodeInt *newLeafNode;
  bufMgr->allocPage(file, newLeafID, newLeaf);
  newLeafNode = (LeafNodeInt *)newLeaf;

  if (leafNode->keyArray[middleIndex] < dataEntry.key && leafOccupancy % 2 == 1){
    middleIndex++;
  }

  for(int i = middleIndex; i < leafOccupancy; i++)  {
	newLeafNode->ridArray[i-middleIndex] = leafNode->ridArray[i];
    newLeafNode->keyArray[i-middleIndex] = leafNode->keyArray[i];
	leafNode->ridArray[i].page_number = INVALID_NUMBER;
    leafNode->keyArray[i] = 0;
  }
  
   if (dataEntry.key < newLeafNode->keyArray[0]){
     childEntry(leafNode, dataEntry);
   }
   else{
     childEntry(newLeafNode, dataEntry);
   }
  pushedUpEntry = new PageKeyPair<int>();
  newLeafNode->rightSibPageNo = leafNode->rightSibPageNo;
  leafNode->rightSibPageNo = newLeafID;

  
  pushedUpEntry->set(newLeafID, newLeafNode->keyArray[0]);
  bufMgr->unPinPage(file, newLeafID, true);
  bufMgr->unPinPage(file, leafNodeID, true);

  if (rootPageNum == leafNodeID) {
    createNewRoot(leafNodeID, pushedUpEntry, true);
  }
}

/*
// The following function is a helper to enter a child for a node. 
//
// @param leafNode, the present value of the leaf where the splitting must occur
// @param insertedData, the data of the child that is trying to be added into the btree
*/
void BTreeIndex::childEntry(LeafNodeInt *leafNode, RIDKeyPair<int> insertedData)
{
  if (leafNode->ridArray[0].page_number == INVALID_NUMBER){
	leafNode->ridArray[0] = insertedData.rid;
    leafNode->keyArray[0] = insertedData.key;
  }
  else{
      int index = 0;
      while(leafNode->ridArray[index].page_number != INVALID_NUMBER) {
        index++;
      }
      while(index >= 1 && (leafNode->keyArray[index-1] > insertedData.key)) {
        leafNode->keyArray[index] = leafNode->keyArray[index-1];
        leafNode->ridArray[index] = leafNode->ridArray[index-1];
        index--;
      }
      leafNode->keyArray[index] = insertedData.key;
      leafNode->ridArray[index] = insertedData.rid;
  }
}

/*
// The follwing method is used to help insert internal nodes into the Btree.
//
// @param nonLeafNode, the node of the non leaf that is being added
// @param insertedInternalNode, the internal node that is tyring to be changed so it can be added
*/
void BTreeIndex::internalNodeEntry(NonLeafNodeInt *nonLeafNode, PageKeyPair<int> *insertedInternalNode)
{ 
  int index = 0;
  while(nonLeafNode->pageNoArray[index+1] != INVALID_NUMBER) {
    index++;
  }

  while(index >= 1 && (nonLeafNode->keyArray[index-1] > insertedInternalNode->key)) {
    nonLeafNode->keyArray[index] = nonLeafNode->keyArray[index-1];
    nonLeafNode->pageNoArray[index+1] = nonLeafNode->pageNoArray[index];
    index--;
  }
  nonLeafNode->keyArray[index] = insertedInternalNode->key;
  nonLeafNode->pageNoArray[index+1] = insertedInternalNode->pageNo;
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

  /**
   * Insert a new entry using the pair <value,rid>.
   * Start from root to recursively find out the leaf to insert the entry in. The insertion may cause splitting of leaf node.
   * This splitting will require addition of new leaf page number entry into the parent non-leaf, which may in-turn get split.
   * This may continue all the way upto the root causing the root to get split. If root gets split, metapage needs to be changed accordingly.
   * Make sure to unpin pages as soon as you can.
   *
   * @param key         Key to insert, pointer to integer/double/char string
   * @param rid         Record ID of a record whose entry is getting inserted into the index.
  **/
void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{
  Page* rootPage;
  PageKeyPair<int> *pushedUpEntry;
  RIDKeyPair<int> insertedData;
  insertedData.set(rid, *((int *)key));
  bufMgr->readPage(file, rootPageNum, rootPage);
  pushedUpEntry = nullptr;
  if(firstRootValue == rootPageNum) {
    insertData(rootPage, rootPageNum, true, insertedData, pushedUpEntry);
  } else {
    insertData(rootPage, rootPageNum, false, insertedData, pushedUpEntry);
  }
}
  /**
   * Helper function to determine if a given key is in a specified range.
   *
   * @param lowVal   Low value of range, pointer to integer / double / char string
   * @param lowOp    Low operator (GT/GTE)
   * @param highVal  High value of range, pointer to integer / double / char string
   * @param highOp   High operator (LT/LTE)
   * @param key      Key to be checked
   * @return true if a given key is in a specified range, false otherwise.
   */
bool BTreeIndex::isKeyFound(int lowVal, const Operator lowOp, int highVal, const Operator highOp, int key)
{
  if(lowOp == GTE && highOp == LTE){
    return key <= highVal && key >= lowVal;
  }
  else if(lowOp == GT && highOp == LTE){
    return key <= highVal && key > lowVal;
  }
  else if(lowOp == GTE && highOp == LT){
    return key < highVal && key >= lowVal;
  }
  else{
    return key < highVal && key > lowVal;
  }
 }

// -----------------------------------------------------------------------------
// BTreeIndex::startScan
// 

  /**
    * Begin a filtered scan of the index.  For instance, if the method is called
    * using ("a",GT,"d",LTE) then we should seek all entries with a value
    * greater than "a" and less than or equal to "d".
    * If another scan is already executing, that needs to be ended here.
    * Set up all the variables for scan. Start from root to find out the leaf page that contains the first RecordID
    * that satisfies the scan parameters. Keep that page pinned in the buffer pool.
    *
    * @param lowVal  Low value of range, pointer to integer / double / char string
    * @param lowOp       Low operator (GT/GTE)
    * @param highVal High value of range, pointer to integer / double / char string
    * @param highOp  High operator (LT/LTE)
    * @throws  BadOpcodesException If lowOp and highOp do not contain one of their their expected values 
    * @throws  BadScanrangeException If lowVal > highval
    * @throws  NoSuchKeyFoundException If there is no key in the B+ tree that satisfies the scan criteria.
  **/
void BTreeIndex::startScan(const void* lowValParm,
           const Operator lowOpParm,
           const void* highValParm,
           const Operator highOpParm)
{
  highValInt = *((int *)highValParm);
  lowValInt = *((int *)lowValParm);

     if((lowOpParm != GTE && lowOpParm != GT) || (highOpParm != LTE && highOpParm != LT)){
    throw BadOpcodesException();
  }
  if(lowValInt > highValInt){
    throw BadScanrangeException();
  }
  highOp = highOpParm;
  lowOp = lowOpParm;

  if(scanExecuting){
    endScan();
  }
  currentPageNum = rootPageNum;
  bufMgr->readPage(file, currentPageNum, currentPageData);
  NonLeafNodeInt* currentScannedNode;
  bool isLeafFound;
  bool isValueFound;
  bool isValidPageFound;
  PageId belowNodeID;
  LeafNodeInt* currentLeafNode;
  int tmpKey;
  if(rootPageNum != firstRootValue)  {
    currentScannedNode = (NonLeafNodeInt *) currentPageData;
    isLeafFound = false;
    while(!isLeafFound){
      currentScannedNode = (NonLeafNodeInt *) currentPageData;
      if(currentScannedNode->level == 1){
        isLeafFound = true;
      }
      belowNodeID = findBelowNode(currentScannedNode, lowValInt);
      bufMgr->unPinPage(file, currentPageNum, false);
      currentPageNum = belowNodeID;
      bufMgr->readPage(file, currentPageNum, currentPageData);
    }
  }
  isValueFound = false;
  while(!isValueFound){
    currentLeafNode = (LeafNodeInt *) currentPageData;
    if(currentLeafNode->ridArray[0].page_number == INVALID_NUMBER) {
      bufMgr->unPinPage(file, currentPageNum, false);
      throw NoSuchKeyFoundException();
    }
    isValidPageFound = true;
  for(int i = 0; i < leafOccupancy && isValidPageFound; i++)
    {
      tmpKey = currentLeafNode->keyArray[i];
      if(currentLeafNode->ridArray[i + 1].page_number == INVALID_NUMBER && i < leafOccupancy - 1){
        isValidPageFound = false;
      }
      
      if(isKeyFound(lowValInt, lowOp, highValInt, highOp, tmpKey)){
        scanExecuting = true;
        isValueFound = true;
		nextEntry = i;
        break;
      }
      else if((highOp == LT && tmpKey >= highValInt) || (highOp == LTE && tmpKey > highValInt)){
        bufMgr->unPinPage(file, currentPageNum, false);
        throw NoSuchKeyFoundException();
      }

    if(!isValidPageFound ||i == leafOccupancy - 1) {
        bufMgr->unPinPage(file, currentPageNum, false);
        if(currentLeafNode->rightSibPageNo == INVALID_NUMBER) {
          throw NoSuchKeyFoundException();
        }
        currentPageNum = currentLeafNode->rightSibPageNo;
        bufMgr->readPage(file, currentPageNum, currentPageData);
      }
    }
  }
}

// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------

/**
  * Fetch the record id of the next index entry that matches the scan.
  * Return the next record from current page being scanned. If current page has been scanned to its entirety, move on to the right sibling of current page, if any exists, to start scanning that page. Make sure to unpin any pages that are no longer required.
  *
  * @param outRid  RecordId of next record found that satisfies the scan criteria returned in this
  * @throws ScanNotInitializedException If no scan has been initialized.
  * @throws IndexScanCompletedException If no more records, satisfying the scan criteria, are left to be scanned.
  **/
void BTreeIndex::scanNext(RecordId& outRid) 
{
  LeafNodeInt* currentLeafNode;
  int tmpKey;
  if(!scanExecuting){
    throw ScanNotInitializedException();
  }
  currentLeafNode = (LeafNodeInt *) currentPageData;
  if(leafOccupancy == nextEntry || (currentLeafNode->ridArray[nextEntry].page_number == INVALID_NUMBER)){
    if(currentLeafNode->rightSibPageNo == INVALID_NUMBER){
      throw IndexScanCompletedException();
    }
    nextEntry = 0;
    bufMgr->unPinPage(file, currentPageNum, false);
    currentPageNum = currentLeafNode->rightSibPageNo;
    bufMgr->readPage(file, currentPageNum, currentPageData);
    currentLeafNode = (LeafNodeInt *) currentPageData;
  }
  tmpKey = currentLeafNode->keyArray[nextEntry];
  if(isKeyFound(lowValInt, lowOp, highValInt, highOp, tmpKey)){
    outRid = currentLeafNode->ridArray[nextEntry];
    nextEntry++;
  }
  else{
    throw IndexScanCompletedException();
  }
}

// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------
//

/**
  * Terminate the current scan. Unpin any pinned pages. Reset scan specific variables.
  *
  * @throws ScanNotInitializedException If no scan has been initialized.
  **/
void BTreeIndex::endScan() {
  if(!scanExecuting) {
    throw ScanNotInitializedException();
  }
  scanExecuting = false;
  bufMgr->unPinPage(file, currentPageNum, false);
  // Reset variable
  nextEntry = -1;
  currentPageNum = static_cast<PageId>(-1);
  currentPageData = nullptr;

}
}
