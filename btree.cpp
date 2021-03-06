/**
 * @author See Contributors.txt for code contributors and overview of BadgerDB.
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
       file = new BlobFile(outIndexName, false);
       headerPageNum = file->getFirstPageNo();
       Page *firstPageInBlob;
       bufMgr->readPage(file, headerPageNum, firstPageInBlob);
       IndexMetaInfo *metaPage = (IndexMetaInfo *)firstPageInBlob;
       rootPageNum = metaPage->rootPageNo;

       // check if index Info matches
       if (relationName != metaPage->relationName || attrByteOffset != metaPage->attrByteOffset
           || attrType != metaPage->attrType)
       {
         throw BadIndexInfoException(outIndexName);
       }

       bufMgr->unPinPage(file, headerPageNum, false);
         //File was not found thus we create a new one
     } catch(FileNotFoundException e) {
       //File did not exist from upon, thus create a new blob file
       file = new BlobFile(outIndexName, true);
       // allocate root and header page
       Page *rootPage;
       Page *firstPage;

       bufMgr->allocPage(file, rootPageNum, rootPage);
       bufMgr->allocPage(file, headerPageNum, firstPage);


       // fill metaPage infor
       IndexMetaInfo *metaPage = (IndexMetaInfo *)firstPage;
       metaPage->rootPageNo = rootPageNum;
       metaPage->attrByteOffset = attrByteOffset;
       metaPage->attrType = attrType;



       relationName.copy(metaPage->relationName, 20, 0);


//       strncpy((char *)(&(metaPage->relationName)), relationName.c_str(), 20);
       //metaPage->relationName[19] = 0;
       firstRootValue = rootPageNum;

       // initiaize root
       LeafNodeInt *root = (LeafNodeInt *)rootPage;
       root->rightSibPageNo = INVALID_NUMBER;

       bufMgr->unPinPage(file, headerPageNum, true);
       bufMgr->unPinPage(file, rootPageNum, true);

       //fill the newly created Blob File using filescan
       FileScan scannedFile(relationName, bufMgr);
       RecordId recordID;
       try
       {
         for(;;)
         {
           scannedFile.scanNext(recordID);
           std::string currentRecord = scannedFile.getRecord();
           insertEntry(currentRecord.c_str() + attrByteOffset, recordID);
         }
       }
       catch(EndOfFileException e)
       {
         // save Btee index file to disk
         bufMgr->flushFile(file);
       }
     }

}

/**
 * BTreeIndex Destructor. 
 * End any initialized scan, flush index file, after unpinning any pinned pages, from the buffer manager
 * and delete file instance thereby closing the index file.
 * Destructor should not throw any exceptions. All exceptions should be caught in here itself. 
 */
BTreeIndex::~BTreeIndex()
{
  scanExecuting = false;
  bufMgr->flushFile(BTreeIndex::file);
  delete file;
  file = nullptr;
}

/**
 * Insert a new entry using the pair <value,rid>. 
 * Start from root to recursively find out the leaf to insert the entry in. The insertion may cause splitting of leaf node.
 * This splitting will require addition of new leaf page number entry into the parent non-leaf, which may in-turn get split.
 * This may continue all the way upto the root causing the root to get split. If root gets split, metapage needs to be changed accordingly.
 * Make sure to unpin pages as soon as you can.
 * @param key     Key to insert, pointer to integer/double/char string
 * @param rid     Record ID of a record whose entry is getting inserted into the index.
**/
const void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{
  RIDKeyPair<int> insertedData;
  insertedData.set(rid, *((int *)key));
  // root
  Page* rootPage;
  // PageId rootPageNum;
  bufMgr->readPage(file, rootPageNum, rootPage);
  PageKeyPair<int> *pushedUpEntry = nullptr;
  if(firstRootValue == rootPageNum) {
    insertData(rootPage, rootPageNum, true, insertedData, pushedUpEntry);
  } else {
    insertData(rootPage, rootPageNum, false, insertedData, pushedUpEntry);
  }
//  insert(rootPage, rootPageNum, firstRootValue == rootPageNum ? true : false, insertedData, pushedUpEntry);
}

/**
 * Helper function to find the next level of page for the key should be in. 
 * @param curPage       The current Page we are checking
 * @param nextNodenum   Return value for the next level page ID
 * @param key           The Key we are checking
*/
//const void BTreeIndex::findNextNonLeafNode(NonLeafNodeInt *nonLeafNode, PageId &belowNodeID, int insertedData)
const PageId BTreeIndex::findBelowNode(NonLeafNodeInt *nonLeafNode, int insertedData)
{
//  int i = nodeOccupancy;
//  while(i >= 0 && (nonLeafNode->pageNoArray[i] == 0))
//  {
//    i--;
//  }
//  while(i > 0 && (nonLeafNode->keyArray[i-1] >= insertedData))
//  {
//    i--;
//  }

  int i = 0;
  if(insertedData < nonLeafNode->keyArray[i]) {
    return nonLeafNode->pageNoArray[i];
  }
  while (i < nodeOccupancy && nonLeafNode->pageNoArray[i+1] != INVALID_NUMBER && (nonLeafNode->keyArray[i] <= insertedData)) {
    i++;
  }

  // while(!(nonLeafNode->keyArray[i] <= insertedData <  nonLeafNode->keyArray[i+1])) {
//      i++;
//  }


  // return nonLeafNode->pageNoArray[i+1];

//  belowNodeID = nonLeafNode->pageNoArray[i];
    return nonLeafNode->pageNoArray[i];
}


/**
 * Recursive function to insert the index entry to the index file
 * @param curPage           The current Page we are checking
 * @param curPageNum        PageId of current Page
 * @param nodeIsLeaf        If the current page is a leaf node or nonleaf node
 * @param dataEntry         Index entry that needs to be inserted
 * @param newchildEntry     A pageKeyPair that contains an entry that is pushed up after splitting a node; it is null if no split in child nodes
*/
const void BTreeIndex::insertData(Page *currentPage, PageId currentNodeID, bool isLeaf, const RIDKeyPair<int> insertedData, PageKeyPair<int> *&pushedUpEntry)
{
  // nonleaf node
  if (isLeaf)
  {
    LeafNodeInt *leafNode = (LeafNodeInt *)currentPage;
    // page is not full
    if (leafNode->ridArray[leafOccupancy - 1].page_number == 0)
    {
      insertLeaf(leafNode, insertedData);
      bufMgr->unPinPage(file, currentNodeID, true);
      pushedUpEntry = nullptr;
    }
    else
    {
      splitLeaf(leafNode, currentNodeID, pushedUpEntry, insertedData);
    }

  }
  else
  {

    NonLeafNodeInt *nonLeafNode = (NonLeafNodeInt *)currentPage;
    // find the right key to traverse
    PageId belowNodeID;
    Page *belowNode;
//    findNextNonLeafNode(nonLeafNode, belowNodeID, insertedData.key);
    belowNodeID = findBelowNode(nonLeafNode, insertedData.key);
    bufMgr->readPage(file, belowNodeID, belowNode);
    // NonLeafNodeInt *nextNode = (NonLeafNodeInt *)nextPage;
    if (nonLeafNode->level == 1) {
      isLeaf = true;
    } else {
      isLeaf = false;
    }
//    isLeaf = nonLeafNode->level == 1;
    insertData(belowNode, belowNodeID, isLeaf, insertedData, pushedUpEntry);
    




    // no split in child, just return
    if (pushedUpEntry == nullptr)
    {
        // unpin current page from call stack
        bufMgr->unPinPage(file, currentNodeID, false);
        return;
    }
//    else
//      {
      // if the curpage is not full
      if (nonLeafNode->pageNoArray[nodeOccupancy] == INVALID_NUMBER)
      {
        // insertData the newchildEntry to curpage
        insertNonLeaf(nonLeafNode, pushedUpEntry);
        pushedUpEntry = nullptr;
        // finish the insert process, unpin current page
        bufMgr->unPinPage(file, currentNodeID, true);
      }
      else
      {
        splitNonLeaf(nonLeafNode, currentNodeID, pushedUpEntry);
      }
//    }
  }
}

/**
 * Recursive function to insert the index entry to the index file
 * @param oldNode           the node that needs to be split
 * @param oldPageNum        PageId of the oldNode
 * @param newchildEntry     A pageKeyPair that contains an entry that is pushed up after splitting a node;
 *                          The value gets updated to contain the new keyPair that needs to be pushed up;
*/
const void BTreeIndex::splitNonLeaf(NonLeafNodeInt *leftNode, PageId leftNodeID, PageKeyPair<int> *&pushedUpEntry)
{
  // allocate a new nonleaf node
  PageId rightNodeID;
  Page *rightPage;
  bufMgr->allocPage(file, rightNodeID, rightPage);
  NonLeafNodeInt *rightNode = (NonLeafNodeInt *)rightPage;

  int middleIndex = nodeOccupancy / 2;
  int pushupIndex = middleIndex;
//  PageKeyPair<int> *pushupEntry = new PageKeyPair<int>();

  // even number of keys
  if (nodeOccupancy % 2 == 0 && leftNode->keyArray[middleIndex] > pushedUpEntry->key) {
    pushupIndex = middleIndex - 1;

  }  else {
   middleIndex++;
   }
//    pushupIndex = pushedUpEntry->key < leftNode->keyArray[middleIndex] ? middleIndex -1 : middleIndex;




//  pushupEntry->set(newPageNum, oldNode>keyArray[pushupIndex]);

//  PageId tmpNewPageNum = rightNodeID;
  int tmpKey = leftNode->keyArray[pushupIndex];

//  middleIndex = pushupIndex + 1;

  // move half the entries to the new node
  for(int i = middleIndex; i < nodeOccupancy; i++)
  {
    rightNode->keyArray[i-middleIndex] = leftNode->keyArray[i];
    rightNode->pageNoArray[i-middleIndex] = leftNode->pageNoArray[i+1];
    leftNode->pageNoArray[i+1] = INVALID_NUMBER;
    leftNode->keyArray[i+1] = 0;
  }


  // remove the entry that is pushed up from current node
  leftNode->keyArray[pushupIndex] = 0;
  leftNode->pageNoArray[pushupIndex] = INVALID_NUMBER;

  rightNode->level = leftNode->level;
  // insertData the new child entry
  if (rightNode->keyArray[0] > pushedUpEntry->key) {
    insertNonLeaf(leftNode, pushedUpEntry);
  } else {
    insertNonLeaf(rightNode, pushedUpEntry);
  }
//  insertNonLeaf(pushedUpEntry->key <  rightNode->keyArray[0] ? leftNode : rightNode, pushedUpEntry);
  // newchildEntry = new PageKeyPair<int>();
//  newchildEntry = pushupEntry;

//   pushedUpEntry->set(tmpNewPageNum, tmpKey);
  pushedUpEntry->set(rightNodeID, tmpKey);
  bufMgr->unPinPage(file, rightNodeID, true);
  bufMgr->unPinPage(file, leftNodeID, true);


  // if the curNode is the root
  if (rootPageNum == leftNodeID) {
    createNewRoot(leftNodeID, pushedUpEntry, false);
    // createNewRoot(leftNode, pushedUpEntry);
  }
}


/**
 * When the root needs to be split, create a new root node and insert the entry pushed up and update the header page 
 * @param firstPageInRoot   The pageId of the first pointer in the root page
 * @param newchildEntry     The keyPair that is pushed up after splitting
*/
const void BTreeIndex::createNewRoot(PageId formerRootID, PageKeyPair<int> *newRoot, bool isLeaf)
// const void BTreeIndex::createNewRoot(NonLeafNodeInt *formerRoot, PageKeyPair<int> *newRoot, bool isLeaf)
{
  // create a new root 
  PageId newRootNodeID;
  Page *newRootPage;
  bufMgr->allocPage(file, newRootNodeID, newRootPage);
  NonLeafNodeInt *newRootNode = (NonLeafNodeInt *)newRootPage;

  // update metadata

  //  newRootNode->level = firstRootValue == rootPageNum ? 1 : 0;
  newRootNode->pageNoArray[0] = formerRootID;
  newRootNode->pageNoArray[1] = newRoot->pageNo;
  newRootNode->keyArray[0] = newRoot->key;
  if (isLeaf) {
    newRootNode->level = 1;
  } else {
    newRootNode->level = 0;
  }

  Page *firstPageInBlob;
  bufMgr->readPage(file, headerPageNum, firstPageInBlob);
  IndexMetaInfo *metaPage = (IndexMetaInfo *)firstPageInBlob;
  metaPage->rootPageNo = newRootNodeID;
  rootPageNum = newRootNodeID;
  // unpin unused page
  bufMgr->unPinPage(file, newRootNodeID, true);
  bufMgr->unPinPage(file, headerPageNum, true);

}

/**
 * Helper function to splitLeafNode when the leafNode is full
 * @param leaf          Leaf node that is full
 * @param leafPageNum   The number of page of that leaf
 * @param newchildEntry The PageKeyPair that need to push up
 * @param dataEntry     The data entry that need to be inserted 
*/
const void BTreeIndex::splitLeaf(LeafNodeInt *leafNode, PageId leafNodeID, PageKeyPair<int> *&pushedUpEntry, const RIDKeyPair<int> dataEntry)
{
  // allocate a new leaf page
  PageId newLeafID;
  Page *newLeaf;
  bufMgr->allocPage(file, newLeafID, newLeaf);
  LeafNodeInt *newLeafNode = (LeafNodeInt *)newLeaf;

  int middleIndex = leafOccupancy/2;
  // odd number of keys
  if (leafNode->keyArray[middleIndex] < dataEntry.key && leafOccupancy % 2 == 1)
  {
    middleIndex++;
  }

  // copy half the page to newLeafNode
  for(int i = middleIndex; i < leafOccupancy; i++)
  {
    newLeafNode->keyArray[i-middleIndex] = leafNode->keyArray[i];
    newLeafNode->ridArray[i-middleIndex] = leafNode->ridArray[i];
    leafNode->keyArray[i] = 0;
    leafNode->ridArray[i].page_number = INVALID_NUMBER;
  }
  
//  if (dataEntry.key > leafNode->keyArray[middleIndex-1])
//  {
//    insertLeaf(newLeafNode, dataEntry);
//  }
//  else
//  {
//    insertLeaf(leafNode, dataEntry);
//  }
   if (dataEntry.key < newLeafNode->keyArray[0])
   {
     insertLeaf(leafNode, dataEntry);
   }
   else
   {
     insertLeaf(newLeafNode, dataEntry);
   }



  // update sibling pointer
  newLeafNode->rightSibPageNo = leafNode->rightSibPageNo;
  leafNode->rightSibPageNo = newLeafID;

  // the smallest key from second page as the new child entry
//  pushedUpEntry = new PageKeyPair<int>();
//  PageKeyPair<int> newKeyPair;
//  newKeyPair.set(newLeafID, newLeafNode->keyArray[0]);
//  pushedUpEntry = &newKeyPair;
  pushedUpEntry = new PageKeyPair<int>();
  pushedUpEntry->set(newLeafID, newLeafNode->keyArray[0]);
  bufMgr->unPinPage(file, newLeafID, true);
  bufMgr->unPinPage(file, leafNodeID, true);


  // if curr page is root
  if (rootPageNum == leafNodeID)
  {
    createNewRoot(leafNodeID, pushedUpEntry, true);
  }
}

const void BTreeIndex::insertLeaf(LeafNodeInt *leafNode, RIDKeyPair<int> insertedData)
{
  // empty leafNode page
  if (leafNode->ridArray[0].page_number == INVALID_NUMBER)
  {
    leafNode->keyArray[0] = insertedData.key;
    leafNode->ridArray[0] = insertedData.rid;
  }
  else
  {

      int index = 0;
      // Find first empty page
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

          // shift insertedData
//          while((leafNode->keyArray[index] > insertedData.key) && index >= 0) {
//            leafNode->keyArray[index+1] = leafNode->keyArray[index];
//            leafNode->ridArray[index+1] = leafNode->ridArray[index];
//            index--;
//          }
//          // insert insertedData
//          leafNode->keyArray[index+1] = insertedData.key;
//          leafNode->ridArray[index+1] = insertedData.rid;



//    int index = leafOccupancy - 1;
//    // find the end
//    while((leafNode->ridArray[index].page_number == INVALID_NUMBER) && index >= 0)
//    {
//      index--;
//    }
//    // shift insertedData
//    while((leafNode->keyArray[index] > insertedData.key) && index >= 0) {
//      leafNode->keyArray[index+1] = leafNode->keyArray[index];
//      leafNode->ridArray[index+1] = leafNode->ridArray[index];
//      index--;
//    }
//    // insert insertedData
//    leafNode->keyArray[index+1] = insertedData.key;
//    leafNode->ridArray[index+1] = insertedData.rid;
  }
}

const void BTreeIndex::insertNonLeaf(NonLeafNodeInt *nonLeafNode, PageKeyPair<int> *insertedInternalNode)
{
  
  int index = 0;
  // Find first empty page
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


//  int i = nodeOccupancy;
//  while(i >= 0 && (nonLeafNode->pageNoArray[i] == 0))
//  {
//    i--;
//  }
//  // shift
//  while( i > 0 && (nonLeafNode->keyArray[i-1] > insertedInternalNode->key))
//  {
//    nonLeafNode->keyArray[i] = nonLeafNode->keyArray[i-1];
//    nonLeafNode->pageNoArray[i+1] = nonLeafNode->pageNoArray[i];
//    i--;
//  }
//  // insert
//  nonLeafNode->keyArray[i] = insertedInternalNode->key;
//  nonLeafNode->pageNoArray[i+1] = insertedInternalNode->pageNo;
}

/**
 * Begin a filtered scan of the index.  For instance, if the method is called 
 * using ("a",GT,"d",LTE) then we should seek all entries with a value 
 * greater than "a" and less than or equal to "d".
 * If another scan is already executing, that needs to be ended here.
 * Set up all the variables for scan. Start from root to find out the leaf page that contains the first RecordID
 * that satisfies the scan parameters. Keep that page pinned in the buffer pool.
 * @param lowVal  Low value of range, pointer to integer / double / char string
 * @param lowOp   Low operator (GT/GTE)
 * @param highVal High value of range, pointer to integer / double / char string
 * @param highOp  High operator (LT/LTE)
 * @throws  BadOpcodesException If lowOp and highOp do not contain one of their their expected values 
 * @throws  BadScanrangeException If lowVal > highval
 * @throws  NoSuchKeyFoundException If there is no key in the B+ tree that satisfies the scan criteria.
**/

const void BTreeIndex::startScan(const void* lowValParm,
           const Operator lowOpParm,
           const void* highValParm,
           const Operator highOpParm)
{
  
  lowValInt = *((int *)lowValParm);
  highValInt = *((int *)highValParm);

//  if(!((lowOpParm == GTE or lowOpParm == GT) && (highOpParm == LTE or highOpParm == LT)))
     if((lowOpParm != GTE && lowOpParm != GT) || (highOpParm != LTE && highOpParm != LT))
  {
    throw BadOpcodesException();
  }
  if(lowValInt > highValInt)
  {
    throw BadScanrangeException();
  }

  lowOp = lowOpParm;
  highOp = highOpParm;

  // Scan is already started
  if(scanExecuting)
  {
    endScan();
  }

  currentPageNum = rootPageNum;
  // Start scanning by reading root page into the buffer pool
  bufMgr->readPage(file, currentPageNum, currentPageData);

  // root is not a leaf
  if(rootPageNum != firstRootValue)
  {
    // Cast
    NonLeafNodeInt* currentScannedNode = (NonLeafNodeInt *) currentPageData;
    bool isLeafFound = false;
    while(!isLeafFound)
    {
      // Cast page to node
      currentScannedNode = (NonLeafNodeInt *) currentPageData;
      // Check if this is the level above the leaf, if yes, the next level is the leaf
      if(currentScannedNode->level == 1)
      {
        isLeafFound = true;
      }

      // Find the leaf
      PageId belowNodeID;
//      findNextNonLeafNode(currentScannedNode, nextPageNum, lowValInt);
      belowNodeID = findBelowNode(currentScannedNode, lowValInt);
      // Unpin
      bufMgr->unPinPage(file, currentPageNum, false);
      currentPageNum = belowNodeID;
      // read the nextPage
      bufMgr->readPage(file, currentPageNum, currentPageData);
    }
  }
  // Now the curNode is leaf node try to find the smallest one that satisfy the OP
  bool isValueFound = false;
  while(!isValueFound){
    // Cast page to node
    LeafNodeInt* currentLeafNode = (LeafNodeInt *) currentPageData;
    // Check if the whole page is null
    if(currentLeafNode->ridArray[0].page_number == INVALID_NUMBER) {
      bufMgr->unPinPage(file, currentPageNum, false);
      throw NoSuchKeyFoundException();
    }
    // Search from the left leaf page to the right to find the fit
//    bool isInvalidPageFound = false;
     bool isValidPageFound = true;
     int tmpKey;
//    for(int i = 0; i < leafOccupancy && !isInvalidPageFound; i++)
  for(int i = 0; i < leafOccupancy && isValidPageFound; i++)
    {
      tmpKey = currentLeafNode->keyArray[i];
      // Check if the next one in the tmpKey is not inserted
      if(i < leafOccupancy - 1 && currentLeafNode->ridArray[i + 1].page_number == INVALID_NUMBER)
      {
//        isInvalidPageFound = true;
        isValidPageFound = false;
      }
      
      if(isKeyFound(lowValInt, lowOp, highValInt, highOp, tmpKey))
      {
        // select
        nextEntry = i;
        scanExecuting = true;
        isValueFound = true;
        break;
      }
      else if((highOp == LTE && tmpKey > highValInt) || (highOp == LT && tmpKey >= highValInt))
      {
        bufMgr->unPinPage(file, currentPageNum, false);
        throw NoSuchKeyFoundException();
      }
      
      // Did not find any matching tmpKey in this leaf, go to next leaf
//      if(i == leafOccupancy - 1 || isInvalidPageFound){
    if(i == leafOccupancy - 1 || !isValidPageFound) {
        //unpin page
        bufMgr->unPinPage(file, currentPageNum, false);
        //did not find the matching one in the more right leaf
        if(currentLeafNode->rightSibPageNo == INVALID_NUMBER) {
          throw NoSuchKeyFoundException();
        }
        currentPageNum = currentLeafNode->rightSibPageNo;
        bufMgr->readPage(file, currentPageNum, currentPageData);
      }
    }
  }
}

/**
  * Fetch the record id of the next index entry that matches the scan.
  * Return the next record from current page being scanned. If current page has been scanned to its entirety, move on to the right sibling of current page, if any exists, to start scanning that page. Make sure to unpin any pages that are no longer required.
  * @param outRid RecordId of next record found that satisfies the scan criteria returned in this
  * @throws ScanNotInitializedException If no scan has been initialized.
  * @throws IndexScanCompletedException If no more records, satisfying the scan criteria, are left to be scanned.
**/
const void BTreeIndex::scanNext(RecordId& outRid) 
{
  if(!scanExecuting)
  {
    throw ScanNotInitializedException();
  }
    // Cast page to node
  LeafNodeInt* currentLeafNode = (LeafNodeInt *) currentPageData;
  if(leafOccupancy == nextEntry || currentLeafNode->ridArray[nextEntry].page_number == INVALID_NUMBER)
  {
    // Unpin page and read papge
//    bufMgr->unPinPage(file, currentPageNum, false);
    // No more next leaf
    if(currentLeafNode->rightSibPageNo == INVALID_NUMBER)
    {
      throw IndexScanCompletedException();
    }
    bufMgr->unPinPage(file, currentPageNum, false);
    currentPageNum = currentLeafNode->rightSibPageNo;
    bufMgr->readPage(file, currentPageNum, currentPageData);
    currentLeafNode = (LeafNodeInt *) currentPageData;
    // Reset nextEntry
    nextEntry = 0;
  }
 
  // Check  if rid satisfy
  int tmpKey = currentLeafNode->keyArray[nextEntry];
  if(isKeyFound(lowValInt, lowOp, highValInt, highOp, tmpKey))
  {
    outRid = currentLeafNode->ridArray[nextEntry];
    // Incrment nextEntry
    nextEntry++;
    // If current page has been scanned to its entirety
  }
  else
  {
    throw IndexScanCompletedException();
  }
}

/**
  * Terminate the current scan. Unpin any pinned pages. Reset scan specific variables.
  * @throws ScanNotInitializedException If no scan has been initialized.
**/
const void BTreeIndex::endScan() {
  if(!scanExecuting) {
    throw ScanNotInitializedException();
  }
  scanExecuting = false;
  // Unpin page
  bufMgr->unPinPage(file, currentPageNum, false);
  // Reset variable
  nextEntry = -1;
  currentPageNum = static_cast<PageId>(-1);
  currentPageData = nullptr;

}

/**
  * Helper function to check if the key is satisfies
  * @param lowVal   Low value of range, pointer to integer / double / char string
  * @param lowOp    Low operator (GT/GTE)
  * @param highVal  High value of range, pointer to integer / double / char string
  * @param highOp   High operator (LT/LTE)
  * @param val      Value of the key
  * @return True if satisfies False if not
  *
**/
// isKeyFound
const bool BTreeIndex::isKeyFound(int lowVal, const Operator lowOp, int highVal, const Operator highOp, int key)
{
  if(lowOp == GTE && highOp == LTE)
  {
    return key <= highVal && key >= lowVal;
  }
  else if(lowOp == GT && highOp == LTE)
  {
    return key <= highVal && key > lowVal;
  }
  else if(lowOp == GTE && highOp == LT)
  {
    return key < highVal && key >= lowVal;
  }
  else
  {
    return key < highVal && key > lowVal;
  }
}

}
