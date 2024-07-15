/*
 *
 * This code is is copied from R and adapted by Peter Langfelder
 *
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1997--2020  The R Core Team
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

#include <Rcpp.h>
#include <R.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace std;

#define NIL -1

/* We use unions here because Solaris gcc -O2 has trouble with
   casting + incrementing pointers.  We use tests here, but R currently
   assumes int is 4 bytes and double is 8 bytes.
 */

union foo {
    double d;
    unsigned int u[2];
};

/*
   Integer keys are hashed via a random number generator
   based on Knuth's recommendations.  The high order K bits
   are used as the hash code.

   NB: lots of this code relies on M being a power of two and
   on silent integer overflow mod 2^32.

   Currently the hash table is implemented as a (signed) integer
    array.  So there are two 31-bit restrictions, the length of the
    array and the values.  The values are initially NIL (-1).  O-based
    indices are inserted by isDuplicated, and invalidated by setting
    to NA_INTEGER.
*/

/*
  Choose M to be the smallest power of 2
  not less than 2*n and set K = log2(M).
  Need K >= 1 and hence M >= 2, and 2^M < 2^31-1, hence n <= 2^29.

  Dec 2004: modified from 4*n to 2*n, since in the worst case we have
  a 50% full table, and that is still rather efficient -- see
  R. Sedgewick (1998) Algorithms in C++ 3rd edition p.606.
*/

template <typename T>
class HashData_Base
{
  protected:
    int K, nomatch;
    size_t M;
    int nmax;
    vector <int> hashTable;
    const T * table;
    int initialized;
    size_t nunique;
    vector <int> unique;

  public:

    size_t scatter(const int key) const {return 3141592653U * key >> (32 - K); }
    virtual int equal(const int , const T ) const = 0;
    virtual size_t hash(const T ) const = 0;

    void MKsetup(int n, const int nmax)
    {
      if(n < 0) /* protect against overflow to -ve */
          stop("length is too large for hashing");
      if (nmax != NA_INTEGER && nmax != 1) n = nmax;
      size_t n2 = 2U * (size_t) n;
      M = 2;
      K = 1;
      while (M < n2) {
          M *= 2;
          K++;
      }
      this->nmax = n;
    }

    void init(const int n, const int nmax, const int nomatch, const T * table)
    {
      MKsetup(n, nmax);
      hashTable.reserve(M);
      fill_n(hashTable.begin(), M, NIL);
      this->nomatch = nomatch;
      this->table = table;
      this->nunique = 0;
      initialized = 1;
    }

    HashData_Base(const int n, const int nmax, const int nomatch, const T * table)
    {
      init(n, nmax, nomatch, table);
    }

    HashData_Base()
    {
      initialized = 0;
    }

    int isDuplicated(const T * x, const size_t indx);
    void removeEntry(const T x);
    void DoHashing(const T * table, const size_t n, const int saveUnique);
    void UndoHashing(const T * x, const size_t n);

    size_t nUnique() { return this->nunique; }
    int isUnique(size_t i) { return this->unique[i]; }
    vector <int> isUnique() { return this->unique; }
};


/* Open address hashing 
   Collision resolution is by linear probing 
   The table is guaranteed large so this is sufficient */
/* Here table is the array of values referenced bythe hash table; x is the value being checked/referenced and indx is the
 * reference to be stored if x is not in the table yet 
*/

template <typename T>
int HashData_Base <T> :: isDuplicated(const T * x, const size_t indx)
{
    if (!initialized)
      stop("HashData_Base instance is not intialized.");

    size_t i = hash(x[indx]);
    while (hashTable[i] != NIL)
    {
      if (equal(hashTable[i], x[indx])) return hashTable[i]>=0 ? 1 : 0;
         // Note to self: hashTable can be negative wither because it is unoccupied (NIL) or because the relevant entry in
         // table is not supposed to be hashed (NA_INTEGER)
      i = (i+1) % M;
    }
    if (nmax-- < 0) stop("hash table is full");
    hashTable[i] = indx;
    return 0;
}

template <typename T>
void HashData_Base <T> :: removeEntry(const T x)
{
    if (!initialized)
      stop("HashData_Base instance is not intialized.");
    size_t i = hash(x);
    while (hashTable[i] >= 0) {
        if (equal(hashTable[i], x)) 
        {
            hashTable[i] = NA_INTEGER;  /* < 0, only index values are inserted */
            return;
        }
        i = (i + 1) % M;
    }
}

/* Build a hash table and optionally save information on duplication */
template <typename T>
void HashData_Base <T> :: DoHashing(const T * table, const size_t n, const int saveUnique)
{
    if (!initialized)
      stop("HashData_Base instance is not intialized.");
    this->nunique = 0;
    if (saveUnique) {
      this->unique.reserve(n);
      fill_n(this->unique.begin(), n, 0);
    } else this->unique.clear();
    for (size_t i = 0; i < n; i++) if (!isDuplicated(table, i)) 
    {
       this->nunique++; 
       if (saveUnique) this->unique[i] = 1;
    }
}

/* invalidate entries: normally few */

template <typename T>
void HashData_Base <T> :: UndoHashing(const T * x, const size_t n)
{
    if (!initialized)
      stop("HashData_Base instance is not intialized.");
    for (size_t i = 0; i < n; i++) removeEntry(x[i]);
}

/* ======================================================================================
 *
 * Specialization for IntegerVector
 *
 * ======================================================================================*/

class HashData_Int: public HashData_Base <int>
{
  public: 
    int equal(const int i, const int x) const
    {  
      if (i < 0) return 0;
      return (this->table[i] == x);
    }

    size_t hash(const int val) const
    {
      if (val==NA_INTEGER) return 0;
      return this->scatter(val);
    }

    IntegerVector HashLookup(IntegerVector x)
    // Note: this returns C indices (0-based). For R need to add 1.
    {
      if (!this->initialized)
        stop("HashData_Base instance is not intialized.");
      size_t n = x.length();
      IntegerVector res(n, this->nomatch);
      //Rprintf("HashLookup: nomatch=%d\n", this->nomatch);
      for (size_t i=0; i<n; i++)
      {
        int ind = hash(x[i]);
        while (this->hashTable[ind] != NIL)
        {
          if (equal(hashTable[ind], x[i]))
                 res[i] = this->hashTable[ind] >= 0 ? this->hashTable[ind] : this->nomatch;
          ind = (ind + 1) % this->M;
        }
        // res[i] = this->nomatch;
      }
      return res;
    }

    size_t intersectSize(IntegerVector x)
    {
      if (!this->initialized)
        stop("HashData_Base instance is not intialized.");
      size_t n = x.length();
      size_t res = 0;
      for (size_t i=0; i<n; i++)
      {
        int ind = hash(x[i]);
        while (this->hashTable[ind] != NIL)
        {
          if (equal(hashTable[ind], x[i]) && this->hashTable[ind] >= 0) { res++; break; };
          ind = (ind + 1) % this->M;
        }
      }
      return res;
    }

    IntegerVector intersect(IntegerVector x)
    {
      if (!this->initialized)
        stop("HashData_Base instance is not intialized.");
      size_t n = x.length(), nInter = 0;
      IntegerVector res( intersectSize(x) );
      
      for (size_t i=0; i<n; i++)
      {
        int ind = hash(x[i]);
        while (this->hashTable[ind] != NIL)
        {
          if (equal(hashTable[ind], x[i]) && hashTable[ind] >= 0)
                 { res[nInter] = this->table[ hashTable[ind] ]; nInter++; break; }
          ind = (ind + 1) % this->M;
        }
      }
      return res;
    }
      
    void DoHashing(const IntegerVector table, const int saveUnique = 0) 
    { 
       HashData_Base <int>::DoHashing ( (int *) table.begin(), size_t(table.length()), saveUnique); 
    }

    void UndoHashing(IntegerVector x) { HashData_Base <int>::UndoHashing(x.begin(), x.length()); }

    void init(const IntegerVector table, const int nomatch, IntegerVector incomparables, int saveUnique = 0)
    {
       this->HashData_Base <int> ::init(size_t(table.length()), 1, nomatch, (int *) table.begin());
       DoHashing(table, saveUnique);
       UndoHashing(incomparables);
    }

    HashData_Int(const IntegerVector table, const int nomatch, IntegerVector incomparables, int saveUnique = 0): 
        HashData_Base(size_t(table.length()), 1, nomatch, (int *) table.begin())
    {
      DoHashing(table, saveUnique);
      UndoHashing(incomparables);
    }

    HashData_Int(): HashData_Base() {};
};

  
/* ======================================================================================
 *
 * Specialization for NumericVector
 *
 * ======================================================================================*/


class HashDataDouble: virtual public HashData_Base <double>
{
  public:
    int equal(const int i, const double x) const
    {  
      if (i < 0) return 0;
      double ti = this->table[i];
      if (!ISNAN(ti) && !ISNAN(x)) return (ti == x);
      else if (ISNA(ti) && ISNA(x)) return 1;
      else if (ISNAN(ti) && ISNAN(x)) return 1;
      else return 0;
    }

    size_t hash(const double val) const
    {
      double tmp = (val == 0.0) ? 0.0 : val;
      if (ISNA(tmp)) tmp = NA_REAL; else if (ISNAN(tmp)) tmp = R_NaN; 
#if 2*SIZEOF_INT == SIZEOF_DOUBLE
      union foo tmpu;
      tmpu.d = tmp;
      return scatter(tmpu.u[0] + tmpu.u[1]);
#else
      return scatter(*((unsigned int *) (&tmp)));
#endif
    }

    NumericVector HashLookup(NumericVector x)
    {
      size_t n = x.length();
      NumericVector res(n);
      for (size_t i=0; i<n; i++)
      {
        int ind = hash(x[i]);
        while (this->hashTable[ind] != NIL)
        {
          if (equal(ind, x[i]))
                 res[i] = this->hashTable[ind] >= 0 ? hashTable[ind] : this->nomatch;
          ind = (ind + 1) % this->M;
        }
        res[i] = this->nomatch;
      }
      return res;
    }
};



