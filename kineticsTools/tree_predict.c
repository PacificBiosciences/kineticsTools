// Description:
// Fast GBM tree predict routintes
//
//  Author: Patrick Marks (pmarks@pacificbiosciences.com)
//

#include <assert.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

int init_native(int a) { return 1; }

#ifndef isnan
inline bool isnan(double x) { return x != x; }
#endif

#ifndef min
inline int min(int x, int y) { return (x > y) ? y : x; }
#endif

/// <summary>
/// Walk the tree for each example, and sum up the leaf nodes.  Emit the total
/// scores for each observation.
/// </summary>
void innerPredict(float radPredF[], float** dataMatrix, int nCtxs, int left[], int right[],
                  int missing[], float splitCode[], int splitVar[], int cSplits[], int varTypes[],
                  float initialValue, int treeSize, int numTrees, int maxCSplitSize)
{

    int tStep = 50;
    int obsStep = 60;

    for (int i = 0; i < nCtxs; i++) {
        radPredF[i] = initialValue;
    }

    for (int t0 = 0; t0 < numTrees; t0 += tStep) {
        for (int obs0 = 0; obs0 < nCtxs; obs0 += obsStep) {
            for (int t = t0; t < min(t0 + tStep, numTrees); t++) {
                int offset = t * treeSize;

                for (int iObs = obs0; iObs < min(obs0 + obsStep, nCtxs); iObs++) {
                    int iCurrentNode = 0;
                    while (splitVar[offset + iCurrentNode] != -1) {
                        float dX = dataMatrix[splitVar[offset + iCurrentNode]][iObs];
                        // missing?
                        if (isnan(dX)) {
                            iCurrentNode = missing[offset + iCurrentNode];
                        }
                        // continuous?
                        else if (varTypes[splitVar[offset + iCurrentNode]] == 0) {
                            if (dX < splitCode[offset + iCurrentNode]) {
                                iCurrentNode = left[offset + iCurrentNode];
                            } else {
                                iCurrentNode = right[offset + iCurrentNode];
                            }
                        } else  // categorical
                        {
                            int iCatSplitIndicator =
                                cSplits[((int)splitCode[offset + iCurrentNode] * maxCSplitSize) +
                                        (int)dX];
                            if (iCatSplitIndicator == -1) {
                                iCurrentNode = left[offset + iCurrentNode];
                            } else if (iCatSplitIndicator == 1) {
                                iCurrentNode = right[offset + iCurrentNode];
                            } else  // categorical level not present in training
                            {
                                iCurrentNode = missing[offset + iCurrentNode];
                            }
                        }
                    }
                    radPredF[iObs] +=
                        (float)splitCode[offset + iCurrentNode];  // add the prediction
                }
            }  // iObs
        }      // iTree
    }
}

static uint32_t modToCanonicalMap[8] = {0, 1, 2, 3, 0, 1, 1, 1};

/// <summary>
/// Walk the tree for each example, and sum up the leaf nodes.  Emit the total
/// scores for each observation.
/// </summary>
void innerPredictCtx(int ctxSize, float radPredF[], uint64_t contextPack[], int nCtxs, int left[],
                     int right[], int missing[], float splitCode[], int16_t splitVar[],
                     int varTypes[], float initialValue, int treeSize, int numTrees,
                     int maxCSplitSize)
{

    // contextPack contains 24 3-bit numbers in feature order

    uint32_t* uintSplitCode = (uint32_t*)splitCode;

    int tStep = 20;
    int obsStep = 1000;

    for (int i = 0; i < nCtxs; i++) {
        radPredF[i] = initialValue;
    }

    for (int t0 = 0; t0 < numTrees; t0 += tStep) {
        for (int obs0 = 0; obs0 < nCtxs; obs0 += obsStep) {
            for (int t = t0; t < min(t0 + tStep, numTrees); t++) {
                int offset = t * treeSize;

                for (int iObs = obs0; iObs < min(obs0 + obsStep, nCtxs); iObs++) {
                    uint64_t ctx = contextPack[iObs];

                    int currentNode = offset;
                    while (splitVar[currentNode] >= 0) {
                        int currentVar = splitVar[currentNode];
                        int ctxPos = currentVar;

                        // Canonical feature means feature over canonical bases A,C,G,T
                        int isCanonicalFeature = currentVar >= ctxSize;

                        if (isCanonicalFeature) ctxPos = currentVar - ctxSize;

                        // context is packed 4 bits per slot, lower 3 bits are the modified base code
                        uint32_t dX = (ctx >> (4 * ctxPos)) & 0x7;

                        if (isCanonicalFeature)
                            // Need the canonical base -- convert
                            // from the general base back to the canonical base
                            dX = modToCanonicalMap[dX];

                        // split code contains packed indicators for each categorical level
                        uint32_t splitPack = uintSplitCode[currentNode];
                        uint32_t ind = (splitPack >> dX) & 0x1;

                        if (ind == 0) {
                            // Left node comes precomputed with offset
                            currentNode = left[currentNode];
                        } else {
                            // Right node come precomputed with offset
                            currentNode = right[currentNode];
                        }
                    }
                    radPredF[iObs] += (float)splitCode[currentNode];  // add the prediction
                }
            }  // iObs
        }      // iTree
    }
}
