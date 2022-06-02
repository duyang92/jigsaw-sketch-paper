#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <random>
#include <unordered_map>
#include <algorithm>
#include <ctime>
#include "MurmurHash3.h"

using namespace std;

#define KEY_SIZE 13
#define DECAY_PARAM 1.08

#define HASH_SEED 1324534127
//static uint32_t HASH_SEED;

static mt19937 rng(time(0));

//for recording information in the unordered_map
struct CmpFunc {
    bool operator()(const uint8_t *keyA, const uint8_t *keyB) const {
        return memcmp(keyA, keyB, KEY_SIZE) == 0;
    }
};

struct HashFunc {
    unsigned int operator()(const uint8_t* key) const {
        unsigned int hashValue=0;
        MurmurHash3_x86_32(key, KEY_SIZE,0, &hashValue);;
        return hashValue;
    }
};

struct HeapNode {
    uint8_t key[KEY_SIZE]={0};
    uint32_t counter=0;
};

template<uint32_t HEAP_MAX_SIZE>
class MinHeap {
private:
    HeapNode array[HEAP_MAX_SIZE];
//    uint32_t maxVal = 0;//record the max value in minHeap
    uint32_t size = 0;//record the current size of minHeap
    unordered_map<uint8_t*, uint16_t, HashFunc, CmpFunc> hashTable;
public:
    MinHeap() {}

    uint32_t getMinVal() {
        if (size > 0) {
            return array[0].counter;
        } else {
            return -1;
        }
    }
//
//    uint32_t getMaxVal() {
//        return maxVal;
//    }

    uint32_t getSize(){
        return size;
    }

    uint8_t* getKey(uint32_t nodeIndex){
        return array[nodeIndex].key;
    }
    uint32_t getValue(uint32_t nodeIndex){
        return array[nodeIndex].counter;
    }

    bool findKey(uint8_t *key, uint32_t &nodeIndex, uint32_t &counter) {
        auto iter=hashTable.find(key);
        if(iter==hashTable.end()){
            return false;
        }else{
            nodeIndex=iter->second;
            counter=array[nodeIndex].counter;
            return true;
        }
    }

    void heapPush(uint8_t *key, uint32_t counter) {
        if (size < HEAP_MAX_SIZE) {
            memcpy(array[size].key, key, KEY_SIZE);
            array[size].counter = counter;
            hashTable[array[size].key]=size;
            heapRepairUp(size);
            size++;

//            maxVal = max(maxVal, counter);
        } else {//heap is full
            if (array[0].counter < counter) {
                auto iter=hashTable.find(array[0].key);
                hashTable.erase(iter);

                memcpy(array[0].key,key,KEY_SIZE);
                hashTable[array[0].key]=0;
                array[0].counter=counter;
                heapRepairDown( 0);

//                maxVal=max(maxVal,counter);
            }
        }
    }

    void updateHeapNode(uint32_t nodeIndex,uint32_t newCounter) {//only consider the case that new counter is larger than old counter
        array[nodeIndex].counter = newCounter;
        heapRepairDown(nodeIndex);
//        maxVal = max(maxVal, newCounter);
    }

private:
    void heapRepairDown(uint32_t nodeIndex) {
        uint32_t l, r;
        uint32_t smallestNodeIdx;

        while (nodeIndex < (size >> 1)) {
            l = (nodeIndex << 1) + 1;
            r = (nodeIndex + 1) << 1;

            smallestNodeIdx = nodeIndex;
            if (l < size && array[l].counter < array[smallestNodeIdx].counter) {
                smallestNodeIdx = l;
            }
            if (r < size && array[r].counter < array[smallestNodeIdx].counter) {
                smallestNodeIdx = r;
            }

            if (smallestNodeIdx != nodeIndex) {
                auto iter=hashTable.find(array[nodeIndex].key);
                hashTable.erase(iter);
                iter=hashTable.find(array[smallestNodeIdx].key);
                hashTable.erase(iter);

                HeapNode temp=array[nodeIndex];
                array[nodeIndex]=array[smallestNodeIdx];
                array[smallestNodeIdx]=temp;

                hashTable[array[nodeIndex].key]=nodeIndex;
                hashTable[array[smallestNodeIdx].key]=smallestNodeIdx;

                nodeIndex=smallestNodeIdx;
            } else {
                break;
            }
        }
    }

    void heapRepairUp(uint32_t nodeIndex) {
        uint32_t parentIdx;
        while (nodeIndex > 0) {
            parentIdx = (nodeIndex - 1) >> 1;

            if (array[nodeIndex].counter < array[parentIdx].counter) {
                auto iter=hashTable.find(array[nodeIndex].key);
                hashTable.erase(iter);
                iter=hashTable.find(array[parentIdx].key);
                hashTable.erase(iter);

                HeapNode temp=array[nodeIndex];
                array[nodeIndex]=array[parentIdx];
                array[parentIdx]=temp;
                hashTable[array[nodeIndex].key]=nodeIndex;
                hashTable[array[parentIdx].key]=parentIdx;

                nodeIndex=parentIdx;
            } else {
                break;
            }
        }
    }
};

struct Bucket {
    uint16_t fingerprint = 0;
    uint32_t counter = 0;
};

//minimum version
template<uint32_t ROW_NUM, uint32_t COL_NUM, uint32_t HEAP_MAX_SIZE>
class HeavyKeeper {
private:
    Bucket bucketArray[ROW_NUM][COL_NUM];
    MinHeap<HEAP_MAX_SIZE> minHeap;
public:
    HeavyKeeper() {};

    void insert(uint8_t *key) {
        uint32_t findedNodeIndex;
        uint32_t counterInMinHeap;
        bool flag;
        flag = minHeap.findKey(key, findedNodeIndex, counterInMinHeap);

        unsigned int hashValue=0;
        MurmurHash3_x86_32(key, KEY_SIZE, HASH_SEED+16516513, &hashValue);

        uint16_t fingerprint = hashValue & 0xffff;

        uint32_t minValInHeap = minHeap.getMinVal();

        Bucket *firstEmptyBucketPtr = nullptr;
        uint32_t counterOfMinBucket = -1;
        Bucket *minBucketPtr = nullptr;

        bool addFlag = false;
        uint32_t toUpdateValue = 0;

        for (uint32_t rowIndex = 0; rowIndex < ROW_NUM; rowIndex++) {
            unsigned int hashValue=0;
            MurmurHash3_x86_32(key, KEY_SIZE, HASH_SEED+9813647*rowIndex, &hashValue);
            uint32_t colIndex = hashValue % COL_NUM;

            Bucket &bucket = bucketArray[rowIndex][colIndex];

            if (bucket.counter == 0) {
                firstEmptyBucketPtr = &bucket;
                break;
            } else if (bucket.fingerprint == fingerprint) {
                if (flag || bucket.counter <= minValInHeap) {
                    bucket.counter++;
                    toUpdateValue = bucket.counter;
                }
                addFlag = true;
                break;
            } else {
                if (bucket.counter < counterOfMinBucket) {
                    counterOfMinBucket = bucket.counter;
                    minBucketPtr = &bucket;
                }
            }
        }

        if (addFlag == false) {
            if (firstEmptyBucketPtr != nullptr) {
                firstEmptyBucketPtr->fingerprint = fingerprint;
                firstEmptyBucketPtr->counter = 1;
                toUpdateValue = 1;
            } else if (minBucketPtr != nullptr) {
                double randomVal = (double)rng()/4294967296;
                if (randomVal < pow(DECAY_PARAM, -1.0 * counterOfMinBucket)) {
                    if (counterOfMinBucket == 1) {
                        minBucketPtr->fingerprint = fingerprint;
                        toUpdateValue = 1;
                    } else {
                        minBucketPtr->counter -= 1;
                    }
                }
            }
        }

        if (flag) {
            if (toUpdateValue > counterInMinHeap) {
                minHeap.updateHeapNode(findedNodeIndex, toUpdateValue);
            }
        } else {
            if (minHeap.getSize() < HEAP_MAX_SIZE or toUpdateValue - minValInHeap == 1) {
                minHeap.heapPush(key, toUpdateValue);
            }
        }
    }
    void getEstimatedFlowSizes(unordered_map<uint8_t *, unsigned int, HashFunc, CmpFunc> &estimatedFlowSizes) {
        for (uint32_t i = 0; i < minHeap.getSize(); i++) {
            uint8_t *tempKey = (uint8_t *) malloc(KEY_SIZE);
            memcpy(tempKey, minHeap.getKey(i), KEY_SIZE);
            estimatedFlowSizes[tempKey] = minHeap.getValue(i);
        }
    }
};


unsigned int ReadInTraces(const char* tracePreFix, uint8_t** keys, unordered_map<uint8_t*, unsigned int, HashFunc, CmpFunc>& actualFlowSizes, unsigned int maxItemNum)
{
    unsigned int count = 0;

    unsigned int countInFile = 0;
    for (int datafileCnt = 0; datafileCnt <= 10; ++datafileCnt){
        countInFile = 0;
        char traceFilePath[100];
        sprintf(traceFilePath, "%s%d.dat", tracePreFix, datafileCnt);
        printf("Start reading %s\n", traceFilePath);

        FILE* fin = fopen(traceFilePath, "rb");
        char temp[KEY_SIZE]{ 0 };
        uint8_t* key;
        while (fread(temp, 1, KEY_SIZE, fin) == KEY_SIZE) {
            key = (uint8_t*)malloc(KEY_SIZE);
            memcpy(key, temp, KEY_SIZE);
            keys[count] = key;
            if (actualFlowSizes.find(key) == actualFlowSizes.end()) {
                actualFlowSizes[key] = 1;
            }
            else {
                actualFlowSizes[key] += 1;
            }
            count++;
            if(count>=maxItemNum){
                printf("The dataset has more than %d items, set a larger value for maxItemNum", maxItemNum);
                exit(-1);
            }

            countInFile++;
            if (countInFile % 5000000 == 0) {
                printf("\thave read %u items in %s, the dataset now has %u items\n", countInFile,traceFilePath, count);
            }
        }
        fclose(fin);
        printf("Finish reading %s (%u items), the dataset now has %u items\n", traceFilePath,countInFile,count);

    }
	return count;
}

bool cmpPairFunc(pair<uint8_t *, int> p1, pair<uint8_t *, int> p2) {
    return p1.second > p2.second;
}

int main()
{
    //prepare dataset
    cout << "prepare dataset" << endl;
    unsigned int maxItemNum = 40*1000000;//max number of items
    uint8_t** keys = (uint8_t**)calloc(maxItemNum, sizeof(uint8_t*));//to store keys of all items
    unordered_map<uint8_t*, unsigned int, HashFunc, CmpFunc> actualFlowSizes;//count the actual flow sizes
    unsigned int actualItemNum = ReadInTraces(R"(../../data/)", keys, actualFlowSizes,maxItemNum);
    cout << "number of items: " << actualItemNum << endl;
    cout << "number of flows: " << actualFlowSizes.size() << endl;
    cout << "*********************" << endl;

    //prepare algorithm
    cout << "prepare algorithm" << endl;

    //parameters
    const uint32_t rowNum = 2;
    //top-1000
    const uint32_t heapMaxSize = 1000;
//    const uint32_t colNum = 3053;//50KB for top-1000
//    const uint32_t colNum = 9077;//100KB for top-1000
//     const uint32_t colNum = 15100;//150KB for top-1000
    const uint32_t colNum = 21124;//200KB for top-1000
//    const uint32_t colNum = 27148;//250KB for top-1000

    //200KB for top2000~5000
//    const uint32_t heapMaxSize = 2000;
//    const uint32_t colNum = 18153;
//    const uint32_t heapMaxSize = 3000;
//    const uint32_t colNum = 15183;
//    const uint32_t heapMaxSize = 4000;
//    const uint32_t colNum = 12212;
//    const uint32_t heapMaxSize = 5000;
//    const uint32_t colNum = 9242;

    HeavyKeeper<rowNum, colNum, heapMaxSize> heavyKeeper;

    //output memory info
    double bucketMem = (2 + 18 / 8.0) * rowNum * colNum / 1024.0;
    //Though we use 32-bit counters, 18 bits is large enough to record the largest flow, thus we use 18 bits as the counter size
    //in 64-bits program, each pointer is 8 bytes, and the key slot is 13 bytes
    double minHeapMem = heapMaxSize * (KEY_SIZE + 18 / 8.0) / 1024.0;//122 bits for each slot
    double hashTableMem = heapMaxSize * (8 + 2) / 1024.0;//80 bits for each item: 64 bit for a pointer, and 16 bit for the slot index of min-heap
    cout << "bucketMem: " << bucketMem << "KB" << endl;
    cout << "minHeapMem: " << minHeapMem << "KB" << endl;
    cout << "hashTableMem: " << hashTableMem << "KB" << endl;
    cout << "totalMem:" << bucketMem + minHeapMem + hashTableMem << "KB" << endl;
    cout << "*********************" << endl;

    //insert items
    cout << "insert items" << endl;
    clock_t time1 = clock();
    for (unsigned int i = 0; i < actualItemNum; i++) {
        heavyKeeper.insert(keys[i]);
    }
    clock_t time2 = clock();

    //calculate throughput
    double numOfSeconds = (double) (time2 - time1) / CLOCKS_PER_SEC;//the seconds using to insert items
    double throughput = (actualItemNum / 1000000.0) / numOfSeconds;
    cout << "use " << numOfSeconds << " seconds" << endl;
    cout << "throughput: " << throughput << " Mpps" << ", each insert operation uses " << 1000.0 / throughput << " ns"
         << endl;
    cout << "*********************" << endl;

    //calculate precision, ARE, AAE
    //get sorted acutal flow sizes and sorted estimated flow sizes
    vector<pair<uint8_t *, unsigned int>> actualFlowSizesVector;
    for (auto iter = actualFlowSizes.begin(); iter != actualFlowSizes.end(); iter++) {
        actualFlowSizesVector.push_back(make_pair(iter->first, iter->second));
    }
    sort(actualFlowSizesVector.begin(), actualFlowSizesVector.end(), cmpPairFunc);

    vector<pair<uint8_t *, unsigned int>> estimatedFlowSizesVector;
    unordered_map<uint8_t *, unsigned int, HashFunc, CmpFunc> estimatedFlowSizes;
    heavyKeeper.getEstimatedFlowSizes(estimatedFlowSizes);

    for (auto iter = estimatedFlowSizes.begin(); iter != estimatedFlowSizes.end(); iter++) {
        estimatedFlowSizesVector.push_back(make_pair(iter->first, iter->second));
    }

    cout << "get top-" << estimatedFlowSizes.size() << " flows" << endl;
    sort(estimatedFlowSizesVector.begin(), estimatedFlowSizesVector.end(), cmpPairFunc);
    unsigned int k = min((unsigned int) heapMaxSize, (unsigned int) estimatedFlowSizes.size());

    //get acutal top-k flows
    unordered_map<uint8_t *, unsigned int, HashFunc, CmpFunc> actualTopKFlowSizes;
    for (unsigned int i = 0; i < k; i++) {
        unsigned int actualSize = actualFlowSizesVector[i].second;
        actualTopKFlowSizes[actualFlowSizesVector[i].first] = actualSize;
    }

    cout << "Top-" << k << " flow detection Metrics" << endl;
    k = min((unsigned int) estimatedFlowSizes.size(), k);
    //top-k precision
    unsigned int TPNumOfTopK = 0;//True positive num
    for (int i = 0; i < k; i++) {
        uint8_t *key = estimatedFlowSizesVector[i].first;
        if (actualTopKFlowSizes.find(key) != actualTopKFlowSizes.end()) {
            TPNumOfTopK++;
        }
    }
    cout << "Precision: " << (double) TPNumOfTopK / k << endl;

    //ARE,AAE of the reported top-k flows
    double totalREOfTopK = 0;//total relative error of top-k flows
    double totalAEOfTopK = 0;
    for (int i = 0; i < k; i++) {
        uint8_t *key = estimatedFlowSizesVector[i].first;
        unsigned int estimatedSize = estimatedFlowSizesVector[i].second;
        unsigned int actualSize = 0;
        if (actualFlowSizes.find(key) != actualFlowSizes.end()) {
            actualSize = actualFlowSizes[key];
        }
        totalAEOfTopK += abs((double) estimatedSize - actualSize);
        double rE = abs((double) estimatedSize - actualSize) / actualSize;
        totalREOfTopK += rE;
    }
    cout << "ARE of reported top-k flows: " << totalREOfTopK / k << endl;
    cout << "AAE of reported top-k flows: " << totalAEOfTopK / k << endl;


    //release resources
    for (auto iter = estimatedFlowSizes.begin(); iter != estimatedFlowSizes.end(); iter++) {
        free(iter->first);
    }
    for (unsigned int i = 0; i < actualItemNum; i++) {
        free(keys[i]);
    }
    free(keys);

    return 0;
}