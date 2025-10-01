#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <sstream>
#include <time.h>
#include <bitset>
#include <vector>
#include <sstream>
#include <fstream>

using namespace std;
typedef unsigned long long int ullint;

//writting the expression in tex format
vector<vector<pair<int,int>>> allGraphs;
vector<ullint> spTr;
vector< vector< pair<ullint,ullint>>> refDen;
vector<int> Sym;
vector< vector< vector<pair<ullint,ullint>>>> matDenNigFinal;
vector<vector<vector<pair<ullint,ullint>>>> matNumSignFinal;
vector<vector<vector<int>>> matCoefFinal;

//simple swap function
void swap(int *arr, int v1, int v2){
    int tps = arr[v1];
    arr[v1] = arr[v2];
    arr[v2] = tps;
}
//Swap nodes without change the order of its positions in L R arrays
void Swap(int *L, int *PosL, int v1, int v2){
    swap(L[v1],L[v2]);
    PosL[L[v1]]=v1;
    PosL[L[v2]]=v2;
    swap(L[PosL[v1]],L[PosL[v2]]);
    swap(PosL[v1],PosL[v2]);
    return;
}
//writing
void WriteUint(const string st, const vector< vector< vector< pair<ullint,ullint>>>> &Mat){
    ofstream myfile(st);

    if (!myfile.is_open()) {
        cerr << "Error: could not open " << st << endl;
        return;
    }

    myfile<<"[";
    for(size_t i = 0; i<Mat.size(); i++){
        if(i!=0)
            myfile<<",[";
        else
            myfile<<"[";
        for(int kk=0; kk<Mat[i].size(); kk++){
            myfile<<"[";
            for(int j=0; j<Mat[i][kk].size(); j++){
                myfile<<"["<<Mat[i][kk][j].first<<","<<Mat[i][kk][j].second<<"]";
                if(j!=Mat[i][kk].size()-1)
                    myfile<<",";
            }
            if(kk==Mat[i].size()-1)
                myfile<<"]";
            else
                myfile<<"],";
        }
        myfile<<"]";
    }
    myfile<<"]";
    myfile.close();
}
//write int vector
void WriteInt(const string st, const vector< vector< vector<int>>> &Mat){
    ofstream myfile(st);
    if (!myfile.is_open()) {
        cerr << "Error: could not open " << st << endl;
        return;
    }

    myfile<<"[";
    for(size_t i = 0; i<Mat.size(); i++){
        if(i!=0)
            myfile<<",[";
        else
            myfile<<"[";
        for(int kk=0; kk<Mat[i].size(); kk++){
            myfile<<"[";
            for(int j=0; j<Mat[i][kk].size(); j++){
                myfile<<Mat[i][kk][j];
                if(j!=Mat[i][kk].size()-1)
                    myfile<<",";
            }
            if(kk==Mat[i].size()-1)
                myfile<<"]";
            else
                myfile<<"],";
        }
        myfile<<"]";
    }
    myfile<<"]";
    myfile.close();
}
//count the number of bits one in a given unsigned integer n
unsigned int countSetBits(ullint n){
    unsigned int cn = 0;
    while(n)
    {
        n^=n&~(n-1);
        cn++;
    }
    return cn;
}
//the sum of all coefficints in cycle
int sumCoefCycle(ullint FC, ullint startEdg, const vector<ullint> &nodes, const vector<ullint> &edgesOut, const vector<pair<int,int>> &edgL, const vector<int> &edgesCoef){
    ullint tail;
    int posCH = __builtin_ctzll(startEdg);
    int nxtNode = edgL[posCH].second, posBR;//next node to be process
    int sum = edgesCoef[posCH];
    FC ^= startEdg;//delete the chorde
    while(FC){
        tail = nodes[nxtNode] & FC;//the next branch
        posBR = __builtin_ctzll(tail);
        if(edgesOut[nxtNode] & tail){//Opposite directions
            sum -= edgesCoef[posBR];
            nxtNode = edgL[posBR].first;
        } else {//Same directions
            sum += edgesCoef[posBR];
            nxtNode = edgL[posBR].second;
        }
        FC ^= tail;
    }
    return sum;
}
//Calculate all the edge coefficients in the graph
void calcEdgCoef(const vector<ullint> &allCycles, const vector<ullint> &edgesOut, const vector<ullint> &nodes, const vector<pair<int,int>> &edgL, vector<int> &edgesCoef){
    int sumCoef, nmberTests = 0;
    bool continueTest = true;
    for(int i = 0; i<2*nodes.size(); i++){
        edgesCoef.push_back(1);
    }
    while(continueTest){
        continueTest = false;
        for(ullint CY: allCycles){
            sumCoef = sumCoefCycle(CY, CY & -CY, nodes, edgesOut, edgL, edgesCoef);
            if(sumCoef == 0){
                edgesCoef[__builtin_ctzll(CY & -CY)] += 2;//incrument the coefficient of an arbitrary edge of the cycle CY by 2 to ensure that edge coefficients remains odds
                continueTest = true;
                //break;//repeat again from zero
            }
        }
        nmberTests++;
        if(nmberTests > 100){
            std::cout<<"Over test "<<std::endl;
            break;
        }
    }
}
//if the cycle contains odd numbers of edges or all the nodes of the cycle have one out edge and one in edge ignore it
bool filtreCycles(ullint CY, const vector<ullint> &edgesOut){
    if(countSetBits(CY) & 1)
        return false;
    ullint tmpOut;
    for(int i = 0; i<edgesOut.size(); i++){
        tmpOut = CY & edgesOut[i];
        if(tmpOut){
            if(tmpOut & (tmpOut-1))
                return true;//this indicate that the node i contains two edges pointed in
        }
    }
    return false;
}
//valiade the generated cycle
bool validCycleTest(ullint CY, const vector<ullint> &nodes){
    vector<ullint> PairEdg;
    ullint doubleEdges;
    for(ullint node: nodes){
        doubleEdges = CY & node;
        if(doubleEdges){
            doubleEdges &= doubleEdges-1;
            if(doubleEdges & (doubleEdges-1))//The node contain more than two edges
                return false;
            else
                PairEdg.push_back(CY & node); //The node contain two pairs of edges
        }
    }
    int nEdgTest = 1;
    bool vesited[PairEdg.size()] = {false};
    ullint index = PairEdg[0];//Initialize to first
    vesited[0] = true;
    while(index)
    {
        for(int i = 1; i<PairEdg.size(); i++){
            if(!vesited[i] && (index&PairEdg[i])){
                index ^= PairEdg[i];
                vesited[i] = true;
                nEdgTest++;
            }
        }
    }
    return nEdgTest == PairEdg.size();
}
//All Combining Possibles
void allCombination(const vector<ullint>& fundamental_cycles, const vector<ullint> &nodes, const vector<ullint> &edgesOut, vector<ullint>& all_Cycles, vector<ullint>& tot_Cycles){
    ullint count = 1ULL << fundamental_cycles.size();
    ullint bitXorCy = 0; //cycle defined by its edges
    int sz = 0;
    for(ullint i = 1; i < count; i++){
        bitXorCy = 0;
        sz = 0;
        for(int j = 0; j < fundamental_cycles.size(); j++){
            if(i & (1ULL << j)){
                sz++;
                bitXorCy ^= fundamental_cycles[j];//Combining cycles
            }
        }
        if(sz>1){
            if(validCycleTest(bitXorCy, nodes)){
                tot_Cycles.push_back(bitXorCy);
                if(filtreCycles(bitXorCy, edgesOut)){
                    all_Cycles.push_back(bitXorCy);
                }
            }
        } else {
            tot_Cycles.push_back(bitXorCy);
            if(filtreCycles(bitXorCy, edgesOut)){
                all_Cycles.push_back(bitXorCy);
            }
        }
    }
}
//Finding the fundamental cycles of undirecetd graph
void findFundamentalCyclesGr(vector<ullint> &fundCycles, ullint spanningTree, ullint coTree, const vector<ullint> &nodes, const vector<pair<int,int>> &edgL){
    ullint tmpGr, branch;
    int pos, k, n = nodes.size();
    while(coTree){
        tmpGr = spanningTree | (coTree & -coTree);//temporal graph
        bool vesited[n] = {false};
        k = 0;
        while(k < n){
            for(int j = 0; j<n; j++){
                if(!vesited[j]){
                    k++;
                    vesited[j] = true;
                    branch = tmpGr & nodes[j];
                    if(branch && ((branch & (branch - 1)) == 0)){//The branch is founded
                        tmpGr ^= branch;//Delete this branch
                        pos = __builtin_ctzll(branch);
                        if(j == edgL[pos].first){//process to the previous vertex
                            j = edgL[pos].second;
                            if(vesited[j]){
                                vesited[j] = false;
                                k--;
                            }
                        }
                        else {
                            j = edgL[pos].first;
                            if(vesited[j]){
                                vesited[j] = false;
                                k--;
                            }
                        }  
                    } else {
                        j++;
                    }
                }
            }
        }
        fundCycles.push_back(tmpGr);//Add this cycle
        coTree &= coTree - 1;
    }
}
//simplify
void simplify(const vector<vector<pair<ullint,ullint>>> &matDenNig, const vector<pair<ullint,ullint>> &matNumSign, vector<vector<pair<ullint,ullint>>> &matDenNigF, vector<vector<pair<ullint,ullint>>> &matNumSignF, vector<vector<int>> &matCoefF, int n){
    int szDD = matDenNig.size();
    // Create a vector of indices to keep track of the original positions
    std::vector<int> positions(szDD);
    for (int i = 0; i < positions.size(); ++i) {
        positions[i] = i;
    }
    // Custom comparator function to sort based on the pairs
    std::sort(positions.begin(), positions.end(), [&](int a, int b) {
        for (size_t i = 0; i < matDenNig[a].size(); ++i) {
            if (matDenNig[a][i].first != matDenNig[b][i].first) {
                return matDenNig[a][i] < matDenNig[b][i];
            }
        }
        return false; // If all elements are equal, return false (no need to swap)
    });
    vector<vector<ullint>> tDen(szDD, vector<ullint>(n-1,0));
    for (int i = 0; i < szDD; i++) {
        for (int j = 0; j < n - 1; j++) {
            tDen[i][j] = matDenNig[positions[i]][j].first;
        }
    }
    matDenNigF.reserve(szDD);
    matNumSignF.resize(szDD);
    matCoefF.resize(szDD);
    int prv = 0, sg = 1, r = 0, pos;
    matDenNigF.push_back(matDenNig[positions[0]]);
    matNumSignF[0].push_back(matNumSign[positions[0]]);
    matCoefF[0].push_back(1);
    for(int i = 1; i<szDD; i++){
        pos = positions[i];
        sg = 1;
        if(tDen[i] != tDen[i-1]){
            prv++;
            matDenNigF.push_back(matDenNig[pos]);
            r = 0;
        } else {
            for(int j = 0; j<n-1; j++)
                if(matDenNig[pos][j].second != matDenNig[positions[i-r-1]][j].second)
                    sg = -sg;
            r++;
        }
        matNumSignF[prv].push_back(matNumSign[pos]);
        matCoefF[prv].push_back(sg);
    }
    matDenNigF.resize(prv+1);
    matNumSignF.resize(prv+1);
    matCoefF.resize(prv+1);
}
//Numerator
void numerator(ullint &sigNum, ullint coTree, const vector< pair<ullint,ullint>> &nigCyclesCy){
    ullint tail;
    for(const auto & ncy : nigCyclesCy){
        tail = ncy.first & coTree;
        if(tail)
            if(((tail & (tail-1)) == 0) && (ncy.second & tail))
                sigNum |= tail;
    }
}
//Denominator
void denominator(ullint spanningTree, const vector<ullint> edgesOut, const vector<ullint>& nodes, const vector<pair<int,int>>& edgL, 
                vector< pair<ullint,ullint>> &denominatorsNig){
    ullint spTr = spanningTree, tail;
    int pos, firstIndex, secondIndex, j = 0, n = nodes.size();
    while(spanningTree){
        tail = spanningTree & nodes[j];
        if (tail && ((tail & (tail - 1)) == 0)) { // Check if the tail is note zero and is a power of 2
            pos = __builtin_ctzll(tail); // Use __builtin_ctzll for 64-bit unsigned long long integers
            firstIndex = edgL[pos].first;
            secondIndex = edgL[pos].second;
            if (firstIndex == j) {
                denominatorsNig[secondIndex].first ^= denominatorsNig[firstIndex].first;
                denominatorsNig[secondIndex].second ^= denominatorsNig[firstIndex].second;
                j = secondIndex;//go to the previous node
            } else {
                denominatorsNig[firstIndex].first ^= denominatorsNig[secondIndex].first;
                denominatorsNig[firstIndex].second ^= denominatorsNig[secondIndex].second;
                j = firstIndex;//go to the previous node
            }
            spanningTree ^= tail; // Remove the tail
        } else {
            j++;
            if(j == n)
                j = 0;
        }
    }
    sort(denominatorsNig.begin(),denominatorsNig.end());
    for(int i = 1; i< n; i++){
        
        denominatorsNig[i].second &= denominatorsNig[i].first;
        if(denominatorsNig[i].second & spTr)
            denominatorsNig[i].second ^= denominatorsNig[i].first;
        denominatorsNig[i-1] = denominatorsNig[i];
    }
    denominatorsNig.pop_back();//delete the last element
}
//fractions
void fractions(const vector<ullint> &allSpanTrees, const vector<pair<ullint,ullint>> &denominatorsNig0, const vector<ullint> &edgesOut, const vector<ullint> &nodes, 
                const vector<pair<int,int>> &edgL, const vector< pair<ullint,ullint>> &nigCyclesCy, int n){
    vector<vector<pair<ullint,ullint>>> matDenNig;
    vector<pair<ullint,ullint>> matNumSign, denominatorsNig;
    vector<int> matCf;
    matDenNig.reserve(allSpanTrees.size());
    matNumSign.reserve(allSpanTrees.size());
    for(const auto& spanningTree: allSpanTrees){
        ullint coTree = ((1ULL<< (2 * n)) - 1) & (~spanningTree);//the co-Tree
        denominatorsNig = denominatorsNig0;
        denominator(spanningTree, edgesOut, nodes, edgL, denominatorsNig);//generate the denominator
        ullint sigNum = 0;
        numerator(sigNum, coTree, nigCyclesCy);//generate the numerator
        matDenNig.push_back(denominatorsNig);//add the results to a global vector of denominators and their negative parts
        matNumSign.emplace_back(coTree,sigNum);//add the results to a global vector of numerators and their signes
        matCf.push_back(1);
    }
    //Simplify and writting
    refDen.push_back(matDenNig[0]);
    vector< vector<pair<ullint,ullint>>> matDenNigF, matNumSignF;
    vector< vector<int>> matCoefF;
    simplify(matDenNig, matNumSign, matDenNigF, matNumSignF, matCoefF, n);//simplify all fractions
    matDenNigFinal.push_back(matDenNigF);
    matNumSignFinal.push_back(matNumSignF);
    matCoefFinal.push_back(matCoefF);
}
//The cartesian product of the compressed spanning trees
void combineBits(const vector<ullint> &numbers, ullint current, int index, vector<ullint> &allSpanTrees) {
    if (index == numbers.size()){
        allSpanTrees.push_back(current);
        return;
    }
    ullint number = numbers[index];
    while (number) {
        combineBits(numbers, current | (number & -number), index + 1, allSpanTrees);
        number &= number - 1; // Delete the rightmost set bit 1
    }
}
//Generate the spanning trees
void spanningTreesGenerator(vector<ullint> &trs, vector< vector<ullint>> &Edg, int k, vector<ullint> &allSpanTrees){
    if(k==0)
        combineBits(trs, 0, 0, allSpanTrees);
    for(int i=0; i<k; i++){
        if(Edg[k][i]){
            trs.push_back(Edg[k][i]);
            for(int j=0; j<i; j++)
                Edg[i][j] |= Edg[k][j];
            spanningTreesGenerator(trs, Edg, k-1, allSpanTrees);
            trs.pop_back();
            for(int j=0; j<i; j++)
                Edg[i][j] ^= Edg[k][j];
        }
    }
}
//The main function to calculate the fractions and its divided differences forms
void init(int *L, int *R, int n){
    vector<pair<int,int>> graph;
    for(int j=0; j<n; j++){
        graph.push_back(make_pair(1+j,1+L[j]));
        graph.push_back(make_pair(1+j,1+R[j]));
    }
    allGraphs.push_back(graph);
    //Initial steps 
    vector<ullint> nodes(n,0), edgesOut(n,0);
    vector<pair<int,int>> edgL(2*n);//labeled edges
    vector< vector<ullint>> Edg(n, vector<ullint>(n, 0));
    for (int j = 0; j < n; j++) {
        j > L[j] ? Edg[j][L[j]] |= 1ULL << (2 * L[j]) : Edg[L[j]][j] |= 1ULL << (2 * L[j]);
        j > R[j] ? Edg[j][R[j]] |= 1ULL << (2 * R[j] + 1) : Edg[R[j]][j] |= 1ULL << (2 * R[j] + 1);
        edgL[2*L[j]] = make_pair(j, L[j]);  
        edgL[2*R[j]+1] = make_pair(j, R[j]);
        edgesOut[j] |= 1ULL << (2 * j);//edges that pointed out from the node j
        edgesOut[j] |= 1ULL << (2 * j + 1);//edges that pointed out from the node j
        nodes[j] |= ((1ULL << (2 * L[j]))|(1ULL << (2 * R[j] + 1))|edgesOut[j]);//All the 4 edges
    }
    vector<ullint> trs, allSpanTrees;
    spanningTreesGenerator(trs, Edg, n-1, allSpanTrees);//Steps 1: Generate all spanning trees   
    ullint refSpanTree = allSpanTrees[0];//Select an arbitrary spanning tree as reference tree
    ullint coTree = ((1ULL<< (2 * n)) - 1) & (~refSpanTree); // Co-tree is the complement of the spanning tree
    spTr.push_back(refSpanTree);
    vector<ullint> FundamentalCycles, allCycles, totCycles;
    findFundamentalCyclesGr(FundamentalCycles, refSpanTree, coTree, nodes, edgL);//First we compute the fundamental cycles of the diagram
    allCombination(FundamentalCycles, nodes, edgesOut, allCycles, totCycles);//Now we generate all cycles
    vector<int> edgesCoef;
    calcEdgCoef(allCycles, edgesOut, nodes, edgL, edgesCoef);//Calculate all edge coefficients
    //The following steps is to find all cycles that contain a negative circulair edges
    int sm;
    vector< pair<ullint,ullint>> nigCyclesCy;
    ullint tilN = 0, cy;
    for(ullint cyy: totCycles){
        tilN = 0;
        cy = cyy;
        while(cy){
            sm = sumCoefCycle(cyy, cy & -cy, nodes, edgesOut, edgL, edgesCoef);
            if(sm > 0)
                tilN |= (cy & -cy);
            cy &= cy -1;
        }
        if(tilN)
            nigCyclesCy.emplace_back(cyy,tilN);
    }
    vector< pair<ullint,ullint>> denominatorsNig0;
    for(int i = 0; i<n; i++){
        denominatorsNig0.emplace_back(nodes[i],edgesOut[i]);
    }
    fractions(allSpanTrees, denominatorsNig0, edgesOut, nodes, edgL, nigCyclesCy, n);
}
//Helper function for DFS
void DFSUtil(int *jdfs, int v, bool visited[], int *VrtxDFS, int *L, int *R){
    visited[v] = true;
    VrtxDFS[(*jdfs)++] = v;
    if (!visited[L[v]])
        DFSUtil(jdfs, L[v], visited, VrtxDFS, L, R);
    if (!visited[R[v]])
        DFSUtil(jdfs, R[v], visited, VrtxDFS, L, R);
}

void SSwap(int *tpL, int *tpPosL, int *syL, int v1, int v2, int sz){
    for(int j=0; j<sz; j++){
        Swap(tpL, tpPosL, v1+j, v2+j);
        swap(syL,v1+j, v2+j);
    }
}

void Reverse(int *L, int *PosL, int *syL, int *Indx, int left, int right, int bg, int sz){
    while (left < right) {
        SSwap(L, PosL, syL, Indx[left+bg], Indx[right+bg], sz);
        left++;
        right--;
    }
    return;
}

void reverse(int *arr, int left, int right){
    while (left < right) {
        swap(arr, left, right);
        left++;
        right--;
    }
    return;
}

bool nextPermutation(int *L, int *PosL, int * syL, int *Indx, int *index, int first, int last, int bg, int sz){
    // Find the first index from the right where nums[i] < nums[i+1]
    if(first == last)//one element
        return false;
    int i = last - 1;
    while (i >= first && index[i] >= index[i + 1]) {
        i--;
    }
    // If no such index is found, the array is in descending order, and it's the last permutation
    if (i == first-1) {
        reverse(index, first, last);
        Reverse(L, PosL, syL, Indx, first, last, bg, sz);
        return false;
    }

    // Find the smallest number to the right of nums[i] that is greater than nums[i]
    int j = last;
    while (index[j] <= index[i]) {
        j--;
    }
    // Swap nums[i] and nums[j]
    swap(index, i, j);
    SSwap(L, PosL, syL, Indx[i+bg], Indx[j+bg], sz);
    // Reverse the subarray to the right of nums[i]
    reverse(index, i+1, last);
    Reverse(L, PosL, syL, Indx, i+1, last, bg, sz);
    return true;
}
// Function to compute factorial
int factorial(int n) {
    int res = 1;
    for(int i = 2; i<=n; i++)
        res *= i;
    return res;
}
// Function to perform a cyclic left shift on a list of integers
void cyclicPermutation(int *L, int *PosL, int bg, int en){
    int temp = L[bg];
    for (int j = bg; j < en; j++) {
        L[j] = L[j+1];
    }
    L[en] = temp;
}

inline bool Mini(int *L1, int *L2, const int n){
    for(int i=0; i<n; i++)
    {
        if(L1[i]==L2[i]) continue;
        if(L1[i]>L2[i])
            return true;
        else
            break;
    }
    return false;
}

inline void ADD(int *L1, int *L2, int n){
    for(int k=0; k<n; k++)
        L1[k]=L2[k];
}

inline bool Equal(int *L1,int* L2, const int n){
    for(int i=0; i<n; i++)
        if(L1[i]!=L2[i])
            return false;
    return true;
}

//Move the first element to the last position in the cycle FindR
void moveFirstToLast(int *arr, int *pos, int * syL, int first, int last) {
    if (last <= first) {
        return;  // No need to move elements in arrays of size 0 or 1
    }

    int firstElement = arr[first], firstElementc = syL[first];
    for (int i = first; i < last; i++) {
        arr[i] = arr[i + 1];
        pos[arr[i + 1]] = i;//update the position of the element arr[i + 1]
    }
    arr[last] = firstElement;
    pos[firstElement] = last;//update the position of the first element
}

void moveLastToFirst(int *arr, int *pos, int *syL, int first, int last) {
    if (last <= first) {
        return; 
    }
    //Move the last element to the first position
    int lastElement = arr[last], lastElementc = syL[last], i;
    for (i = last; i > first; i--) {
        arr[i] = arr[i - 1];
        //syL[i] = syL[i-1];
        pos[arr[i - 1]] = i;
    }
    arr[first] = lastElement;
    //syL[first] = lastElementc;
    pos[lastElement] = first;
}

void moveLastToFirstV(int *arr, int *pos, int *syL, int first, int last) {
    if (last <= first) {
        return; 
    }
    //Move the last element to the first position
    int lastElement = arr[last], lastElementc = syL[last], i;
    for (i = last; i > first; i--) {
        arr[i] = arr[i - 1];
        syL[i] = syL[i-1];
        pos[arr[i - 1]] = i;
    }
    arr[first] = lastElement;
    syL[first] = lastElementc;
    pos[lastElement] = first;
    //incriment the cycle elements by one
    for (i = first; i <= last; i++) {
        arr[pos[i]]++;
        if(arr[pos[i]] > last)
            arr[pos[i]] = first;
    }
    //update the positions of all allements
    lastElement = pos[last];
    for (i = last; i > first; i--) {
        pos[i] = pos[i - 1];
    }
    pos[first] = lastElement;
}

void reverseV(int *arr, int *pos, int *syL, int first, int last){
    while (first<last){
        swap(arr[first], arr[last]);
        swap(syL[first], syL[last]);
        pos[arr[first]] = first;
        pos[arr[last]] = last;
        swap(arr[pos[first]], arr[pos[last]]);
        swap(pos[first],pos[last]);
        first++;
        last--;
    }
    return;
}
bool firstVisite = true;
void disjCyclesCycGen(int *partition, int len, int n, int eps){
    int Indx[n+1], count[n+1], LisGenDet[n], LisSizDet[n], matbeg[n];
    int L[n], R[n], syL[n], PosL[n], PosR[n], index[n][n], Lev[n], Rev[n], first[n], last[n];
    for (int i = 0; i <= n; i++) {
        Indx[i] = 0;
        count[i] = 0;
        if(i >= len)
            partition[i] = 0;
    }
    // Count occurrences of each integer
    for (int i = 0; i < n; i++) {
        L[i] = i;
        R[i] = i;
        syL[i] = i;
        PosL[i] = i;
        PosR[i] = i;
        count[partition[i]]++;
    }
    for(int r=1; r<=len; r++){
        Indx[r] = partition[r-1];
        Indx[r] += Indx[r-1]; 
    }
    
    int rr = 0, j = 1, idx = 0;
    for (int i = 1; i <= n; i++) {
        if(count[i]>1){
            LisGenDet[rr]=count[i];
            for(int k = 0; k<count[i]; k++){
                index[rr][k] = k;
            }
            LisSizDet[rr] = i;
            matbeg[rr]=j-1;
            rr++;
        }
        j = j+count[i];
        Lev[i-1] = i;
        Rev[i-1] = i;
    }
    Lev[n-1] = 0;
    Rev[n-1] = 0;
    int allDet = rr;
    // Calculate the total number of cyclic permutations
    int totalCombinations = 1;
    for (int i = 0; i < len; i++) {
        totalCombinations *= 2*partition[i];
    }
    int sig = 1;
    if(eps == -1)
        for(int k=0; k<len; k++)
            sig=sig*eps;
    else
        sig = (n % 2 == 0) ? 1 : -1;
    
    int var[allDet], varEnd[allDet], totalPermutations = 1;
    if(allDet != 0){
        
        for(int i = 0; i<allDet; i++){
            var[i] = 0;
            varEnd[i] = factorial(LisGenDet[i]);
            totalPermutations *= varEnd[i];
        }
    }
    
    int shifts[len];
    for (int i = 0; i < len; i++) {
        shifts[i] = 0;
    }
    
    do {
        bool hf, eqDl;
        int jdfs = 0;
        for(int k=0; k<len; k++){//R will be the moving of the end element to the begining position of each cycle in L
            R[Indx[k]]=L[Indx[k+1]-1];
            for(int j = Indx[k]+1; j < Indx[k+1]; j++)
                R[j]=L[j-1];
        }

        if(len != 1)
        {
            bool visited[n];
            for(int k=0; k<n; k++)
                visited[k] = false;
            int VrtxDFS[n];
            DFSUtil(&jdfs, 0, visited, VrtxDFS, L, R);
        }

        if(jdfs == n || len == 1)/*if not disconnected diagrams*/
        {
            hf = false;
            for(int k=0; k<n; k++){
                if(R[k]==k||L[k]==k){
                    hf = true;
                    break;
                }
            }
            if(!hf){
                bool canonic = true;//test if the selected diagram labeled in canonic form or not
                int dv = 0, sym;
                if(Equal(L,Lev,n) && Equal(R,Rev,n)){
                    dv = (1ULL<< n )*n;
                } else {
                    int tpL[n], tpPosL[n];
                    int tpR[n], tpPosR[n];
                    for(int k=0; k<n; k++){
                        PosR[R[k]]=k;
                        syL[k] = k;
                    }
                    for(int k=0; k<n; k++){
                        PosL[L[k]] = k;
                    }
                    ADD(tpL,L,n);
                    ADD(tpPosL,PosL,n);
                    int rot[len] = {0};
                    for(int k = 0; k<len; k++)
                        rot[k] = 0;
                    dv = 0;
                    for (int i = 0; i < len; i++) {
                        shifts[i] = 0;
                    }
                    for(int comb = 0; comb < totalCombinations; comb++){// Iterate through all possible combinations
                        if(!canonic) {
                            break;      
                            }
                        if(allDet !=0 ){
                            for(int k = 0; k<allDet; k++)
                                var[k] = 0;
                            for(int per = 0; per < totalPermutations; per++){// Iterate through all possible permutations
                                if(Mini(L,tpL,n)){
                                    canonic = false;   
                                    break;
                                }
                                else if(Equal(L,tpL,n)){
                                    dv++;  
                                }
                                for(int i = allDet-1; i >= 0; i--){
                                    nextPermutation(tpL, tpPosL, syL, Indx, index[i], 0, LisGenDet[i]-1, matbeg[i], LisSizDet[i]);
                                    var[i]++;
                                    if(var[i] == varEnd[i]){
                                        var[i] = 0;
                                    } else {
                                        break;
                                    }
                                }
                            }
                        }
                        
                        else{
                            if(Mini(L,tpL,n)){
                                canonic = false;
                                break;
                            } 
                            else if(Equal(L,tpL,n)){
                            dv++;  
                            }   
                        }
                        for (int i = 0; i < len; i++){// Apply cyclic rotation to each list
                        if(rot[i])
                            moveLastToFirst(tpL, tpPosL, syL, Indx[i], Indx[i+1]-1);
                            moveLastToFirstV(tpL, tpPosL, syL, Indx[i], Indx[i+1]-1);
                            shifts[i]++;
                            if(shifts[i] == partition[i]){
                                if(!rot[i])
                                    rot[i] = 1;
                                else
                                    rot[i] = 0;
                                reverseV(tpL, tpPosL, syL, Indx[i], Indx[i+1]-1);
                            }
                            
                        if(rot[i])
                            moveFirstToLast(tpL, tpPosL, syL, Indx[i], Indx[i+1]-1);//FindR
                            
                            if (shifts[i] == 2*partition[i]) {
                                shifts[i] = 0;
                                
                            } else {
                                break;
                            }
                        }
                    }
                }
                if(canonic){ //Save the results
                    sym = -sig*dv;
                    Sym.push_back(sym);
                    init(L, R, n);
                }
            }
        }
    }while(next_permutation(L, L+n));
}
//Integer partition
void generatePartitions(int n, int nn, int* partition, int currentIndex, int eps) {
    if (nn == 0) {
            disjCyclesCycGen(partition, currentIndex, n, eps);
        return;
    }
    int start = (currentIndex == 0) ? 1 : partition[currentIndex - 1];
    for (int i = start; i <= nn; i++) {
        partition[currentIndex] = i;
        generatePartitions(n, nn - i, partition, currentIndex + 1, eps);
    }
}
// Main function
int main() {
    int n;
    n = 4;
    cout<<"Enter n";
    cin >> n;
    int partition[n];
    cout << "Enter -1 for fermions and +1 bosons";
    int eps = -1;//-1 for Fermion and +1 for Bosons
    cin >> eps;
    generatePartitions(n, n, partition, 0, eps);
    ostringstream stm ;
    stm << n;
    string sn = stm.str();
    string filename1 = "matDenNigF-"+sn + ".dat";
    string filename2 = "matNumSignF-"+sn + ".dat";
    string filename3 = "matCoefF-"+sn + ".dat";
    WriteUint(filename1, matDenNigFinal);
    WriteUint(filename2, matNumSignFinal);
    WriteInt(filename3, matCoefFinal);

    ofstream SPout("spanningTrees-"+sn+".dat");
    SPout << "[";
    for (size_t g = 0; g < spTr.size(); g++) {
        SPout << spTr[g];
        if (g + 1 < spTr.size()) SPout << ",";
    }
    SPout << "]";
    SPout.close();

    ofstream SYM("symmetries-"+sn+".dat");
    SYM << "[";
    for (size_t g = 0; g < Sym.size(); g++) {
        SYM << Sym[g];
        if (g + 1 < Sym.size()) SYM << ",";
    }
    SYM << "]";
    SYM.close();

    ofstream fout("graphs-"+sn+".dat");
    fout << "[";
    for (size_t g = 0; g < allGraphs.size(); g++) {
        fout << "[";
        for (size_t i = 0; i < allGraphs[g].size(); i++) {
            fout << "(" << allGraphs[g][i].first << "," << allGraphs[g][i].second << ")";
            if (i + 1 < allGraphs[g].size()) fout << ", ";
        }
        fout << "]";
        if (g + 1 < allGraphs.size()) fout << ",";
    }
    fout << "]";
    fout.close();

    ofstream refDEN("refDenominators-"+sn+".dat");
    refDEN << "[";
    for(size_t g = 0; g < refDen.size(); g++) {
        refDEN << "[";
        for (size_t i = 0; i < refDen[g].size(); i++) {
            refDEN << "[" << refDen[g][i].first << "," << refDen[g][i].second << "]";
            if (i + 1 < refDen[g].size()) refDEN << ", ";
        }
        refDEN << "]";
        if (g + 1 < refDen.size()) refDEN << ",";
    }           
    refDEN << "]";
    refDEN.close();
    
    ofstream refN("N.dat");
    refN << n ;
    refN.close();
    return 0;
}
