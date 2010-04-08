#include "DynSuffixArray.h"
#include <iostream>

namespace Moses {

DynSuffixArray::DynSuffixArray() {
  m_SA = new vuint_t();
  m_ISA = new vuint_t();
  m_F = new vuint_t();
  m_L = new vuint_t();
  std::cerr << "DYNAMIC SUFFIX ARRAY CLASS INSTANTIATED" << std::endl;
}

DynSuffixArray::~DynSuffixArray() {
  delete m_SA;
  delete m_ISA;
  delete m_F;
  delete m_L;
}

DynSuffixArray::DynSuffixArray(vuint_t* crp) {
  // make native int array and pass to SA builder
  m_corpus = crp; 
  int size = m_corpus->size();
  int* tmpArr = new int[size];
  for(int i=0 ; i < size; ++i) tmpArr[i] = i; 
  
  qsort(tmpArr, 0, size-1);
  
  m_SA = new vuint_t(tmpArr, tmpArr + size);
  //std::cerr << "printing SA " << std::endl;
  //for(int i=0; i < size; ++i) std::cerr << m_SA->at(i) << std::endl;
  delete[] tmpArr;
  std::cerr << "DYNAMIC SUFFIX ARRAY CLASS INSTANTIATED WITH SIZE " << size << std::endl;
  buildAuxArrays();
  //printAuxArrays();
}

void DynSuffixArray::buildAuxArrays() {
  int size = m_SA->size();
  m_ISA = new vuint_t(size);
  m_F = new vuint_t(size);
  m_L = new vuint_t(size);
  
  for(int i=0; i < size; ++i) {
    m_ISA->at(m_SA->at(i)) = i;
    //(*m_ISA)[(*m_SA)[i]] = i;
    (*m_F)[i] = (*m_corpus)[m_SA->at(i)];
    (*m_L)[i] = (*m_corpus)[(m_SA->at(i) == 0 ? size-1 : m_SA->at(i)-1)]; 
  }
}

int DynSuffixArray::rank(unsigned word, unsigned idx) {
/* use Gerlach's code to make rank faster */
  // the number of word in L[0..i]
  int r(0);
  for(unsigned i=0; i < idx; ++i)
    if(m_L->at(i) == word) ++r;
  return r;
}

/* count function should be implemented
 * with binary search over suffix array!! */
int DynSuffixArray::F_firstIdx(unsigned word) {
  // return index of first row where word is found in m_F
  int low = std::lower_bound(m_F->begin(), m_F->end(), word) - m_F->begin();
  if(m_F->at(low) == word) return low;
  else return -1;
}

/* uses rank() and c() to obtain the LF function */
int DynSuffixArray::LF(unsigned L_idx) {
  int fIdx(-1);
  unsigned word = m_L->at(L_idx);
  if((fIdx = F_firstIdx(word)) != -1) 
    return fIdx + rank(word, L_idx);
  else return -1;
}

void DynSuffixArray::insertFactor(vuint_t* newSent, unsigned newIndex) {
  // for sentences
  //stages 1, 2, 4 stay same from 1char case
  //(use last word of new text in step 2 and save Ltmp until last insert?)
  //stage 3...all words of new sentence are inserted backwards
  // stage 2: k=ISA[newIndex], tmp= L[k], L[k]  = newChar
  assert(newIndex <= m_SA->size());
  int k(-1), kprime(-1);
  k = (newIndex < m_SA->size() ? m_ISA->at(newIndex) : m_ISA->at(0)); // k is now index of the cycle that starts at newindex
  int true_pos = LF(k); // track cycle shift (newIndex - 1)
  int Ltmp = m_L->at(k);
  m_L->at(k) = (*newSent)[newSent->size()-1];  // cycle k now ends with correct word
  
  for(int j = newSent->size()-1; j > -1; --j) {
    kprime = LF(k);  // find cycle that starts with (newindex - 1)
    //kprime += ((m_L[k] == Ltmp) && (k > isa[k]) ? 1 : 0); // yada yada
    // only terminal char can be 0 so add new vocab at end
    kprime = (kprime > 0 ? kprime : m_SA->size());  
    true_pos += (kprime <= true_pos ? 1 : 0); // track changes
    // insert everything
    m_F->insert(m_F->begin() + kprime, (*newSent)[j]);
    int theLWord = (j == 0 ? Ltmp : (*newSent)[j-1]);
  
    m_L->insert(m_L->begin() + kprime, theLWord);
    piterate(m_SA, itr) 
      if(*itr >= newIndex) ++(*itr);
  
    m_SA->insert(m_SA->begin() + kprime, newIndex);
    piterate(m_ISA, itr)
      if((int)*itr >= kprime) ++(*itr);
  
    m_ISA->insert(m_ISA->begin() + newIndex, kprime);
    k = kprime;
  }
  // Begin stage 4
  reorder(true_pos, LF(kprime)); // actual position vs computed position of cycle (newIndex-1)
}

void DynSuffixArray::reorder(unsigned j, unsigned jprime) {
  printf("j=%d\tj'=%d\n", j, jprime);
  while(j != jprime) {
    printf("j=%d\tj'=%d\n", j, jprime);
    int tmp, isaIdx(-1);
    int new_j = LF(j);
    // for SA, L, and F, the element at pos j is moved to j'
    tmp = m_L->at(j); // L
    m_L->at(j) = m_L->at(jprime);
    m_L->at(jprime) = tmp;
    tmp = m_SA->at(j);  // SA
    m_SA->at(j) = m_SA->at(jprime);
    m_SA->at(jprime) = tmp;
  
    // all ISA values between (j...j'] decremented
    for(size_t i = 0; i < m_ISA->size(); ++i) {
      if((m_ISA->at(i) == j) && (isaIdx == -1)) 
        isaIdx = i; // store index of ISA[i] = j
      if((m_ISA->at(i) > j) && (m_ISA->at(i) <= jprime)) --(*m_ISA)[i];
    }
    // replace j with j' in ISA
    //isa[isaIdx] = jprime;
    m_ISA->at(isaIdx) = jprime;
    j = new_j;
    jprime = LF(jprime);
  }
}

void DynSuffixArray::deleteFactor(unsigned index, unsigned num2del) {
  int ltmp = m_L->at(m_ISA->at(index));
  int true_pos = LF(m_ISA->at(index)); // track cycle shift (newIndex - 1)
  for(size_t q = 0; q < num2del; ++q) {
    int row = m_ISA->at(index); // gives the position of index in SA and m_F
    std::cerr << "row = " << row << std::endl;
    std::cerr << "SA[r]/index = " << m_SA->at(row) << "/" << index << std::endl;
    true_pos -= (row <= true_pos ? 1 : 0); // track changes
    m_L->erase(m_L->begin() + row);     
    m_F->erase(m_F->begin() + row);     
  
    m_ISA->erase(m_ISA->begin() + index);  // order is important     
    piterate(m_ISA, itr)
      if((int)*itr > row) --(*itr);
  
    m_SA->erase(m_SA->begin() + row);
    piterate(m_SA, itr)
      if(*itr > index) --(*itr);
  }
  m_L->at(m_ISA->at(index))= ltmp;
  reorder(LF(m_ISA->at(index)), true_pos);
  printAuxArrays();
}

void DynSuffixArray::substituteFactor(vuint_t* newSents, unsigned newIndex) {
  std::cerr << "NEEDS TO IMPLEMENT SUBSITITUTE FACTOR\n";
  return;
}

bool DynSuffixArray::getCorpusIndex(const vuint_t* phrase, vuint_t* indices) {
  pair<vuint_t::iterator,vuint_t::iterator> bounds;
  indices->clear();
  size_t phrasesize = phrase->size();
  // find lower and upper bounds on phrase[0]
  bounds = std::equal_range(m_F->begin(), m_F->end(), phrase->at(0));
  // bounds holds first and (last + 1) index of phrase[0] in m_SA
  size_t lwrBnd = size_t(bounds.first - m_F->begin());
  size_t uprBnd = size_t(bounds.second - m_F->begin());
  if(uprBnd - lwrBnd == 0) return false;  // not found
  if(phrasesize == 1) {
    for(size_t i=lwrBnd; i < uprBnd; ++i) {
      indices->push_back(m_SA->at(i));
    }
    return (indices->size() > 0);
  }
  //find longer phrases if they exist
  for(size_t i = lwrBnd; i < uprBnd; ++i) {
    size_t crpIdx = m_SA->at(i);
    if((crpIdx + phrasesize) >= m_corpus->size()) continue; // past end of corpus
    for(size_t pos = 1; pos < phrasesize; ++pos) { // for all following words
      if(m_corpus->at(crpIdx + pos) != phrase->at(pos)) {  // if word doesn't match
        if(indices->size() > 0) i = uprBnd;  // past the phrases since SA is ordered
        break; 
      }
      else if(pos == phrasesize-1) { // found phrase 
        indices->push_back(crpIdx + pos);  // store rigthmost index of phrase 
      }
    }  
  }
  //cerr << "Total count of phrase = " << indices->size() << endl;
  return (indices->size() > 0);
}

void DynSuffixArray::save(FILE* fout) {
  fWriteVector(fout, *m_SA);
}

void DynSuffixArray::load(FILE* fin) {
  fReadVector(fin, *m_SA);
}

int DynSuffixArray::compare(int pos1, int pos2, int max) {
  for (size_t i = 0; i < (unsigned)max; ++i) {
    if((pos1 + i < m_corpus->size()) && (pos2 + i >= m_corpus->size()))
      return 1;
    if((pos2 + i < m_corpus->size()) && (pos1 + i >= m_corpus->size()))
      return -1;
  
    int diff = m_corpus->at(pos1+i) - m_corpus->at(pos2+i);
    if(diff != 0) return diff;
  }
  return 0;
}

void DynSuffixArray::qsort(int* array, int begin, int end) {
  if(end > begin) 
  {
    int index; 
    {	
      index = begin + (rand() % (end - begin + 1));
      int pivot = array[index];
      {
        int tmp = array[index];
        array[index] = array[end];
        array[end] = tmp;
      }
      for(int i=index=begin; i < end; ++i) {
        if (compare(array[i], pivot, 20) <= 0) {
          {
            int tmp = array[index];
            array[index] = array[i];
            array[i] = tmp;
            index++;
          }
        }
      }
      {
        int tmp = array[index];
        array[index] = array[end];
        array[end] = tmp;
      }
    }
    qsort(array, begin, index - 1);
    qsort(array, index + 1,  end);
  }
}



} // end namespace
