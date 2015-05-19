////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// Context.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Context.h,v 1.12 2008-09-08 15:56:18 fgygi Exp $

#ifndef CONTEXT_H
#define CONTEXT_H

#include <iosfwd>
#include <string>
#include <complex>

#if USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

class Context
{
  private:

  struct ContextImpl* pimpl_;

  public:

  int ictxt() const;
  int myrow() const;
  int mycol() const;
  int nprow() const;
  int npcol() const;

  // number of processes in the context
  // returns -1 if current process is not part of this context
  int size() const;
  // position of current process in row-major order
  // returns -1 if current process is not part of this context
  int myproc() const;
  int mype() const;
  // process id of process (irow,icol)
  int pmap(int irow, int icol) const;

  bool onpe0(void) const;
  bool active(void) const;
  void abort(int ierr) const;
  void barrier(void) const;
  void barrier(char scope) const;

  
 
  // double communications
  void dsend(int m, int n, double* a, int lda, int rdest, int cdest) const;
  void drecv(int m, int n, double* a, int lda, int rsrc, int csrc) const;

  void dsum(char scope, char topology,
            int m, int n, double* a, int lda, int rdest, int cdest) const;
  void dsum(char scope, int m, int n, double* a, int lda) const;
  void dsum(int m, int n, double* a, int lda) const;

  void dmax(char scope, char topology,
            int m, int n, double* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const;
  void dmax(char scope, char topology,
            int m, int n, double* a, int lda, int rdest, int cdest) const;
  void dmax(char scope, int m, int n, double* a, int lda) const;
  void dmax(int m, int n, double* a, int lda) const;

  void dmin(char scope, char topology,
            int m, int n, double* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const;
  void dmin(char scope, char topology,
            int m, int n, double* a, int lda, int rdest, int cdest) const;
  void dmin(char scope, int m, int n, double* a, int lda) const;
  void dmin(int m, int n, double* a, int lda) const;

  void dbcast_send(char scope, char topology,
                   int m, int n, double* a, int lda) const;
  void dbcast_send(char scope, int m, int n, double* a, int lda) const;
  void dbcast_send(int m, int n, double* a, int lda) const;

  void dbcast_recv(char scope, char topology,
               int m, int n, double* a, int lda, int rsrc, int csrc) const;
  void dbcast_recv(char scope, int m, int n, double* a,
                   int lda,int rsrc, int csrc) const;
  void dbcast_recv(int m, int n, double* a, int lda,int rsrc, int csrc) const;

  // integer communications
  void isend(int m, int n, int* a, int lda, int rdest, int cdest) const;
  void irecv(int m, int n, int* a, int lda, int rsrc, int csrc) const;
  void isum(char scope, char topology,
            int m, int n, int* a, int lda, int rdest, int cdest) const;
  void isum(char scope, int m, int n, int* a, int lda) const;
  void isum(int m, int n, int* a, int lda) const;

  void imax(char scope, char topology,
            int m, int n, int* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const;
  void imax(char scope, char topology,
            int m, int n, int* a, int lda, int rdest, int cdest) const;
  void imax(char scope, int m, int n, int* a, int lda) const;
  void imax(int m, int n, int* a, int lda) const;

  void imin(char scope, char topology,
            int m, int n, int* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const;
  void imin(char scope, char topology,
            int m, int n, int* a, int lda, int rdest, int cdest) const;
  void imin(char scope, int m, int n, int* a, int lda) const;
  void imin(int m, int n, int* a, int lda) const;

  void ibcast_send(char scope, char topology,
                   int m, int n, int* a, int lda) const;
  void ibcast_send(char scope, int m, int n, int* a, int lda) const;
  void ibcast_send(int m, int n, int* a, int lda) const;

  void ibcast_recv(char scope, char topology,
                   int m, int n, int* a, int lda, int rsrc, int csrc) const;
  void ibcast_recv(char scope, int m, int n,
                   int* a, int lda,int rsrc, int csrc) const;
  void ibcast_recv(int m, int n, int* a, int lda,int rsrc, int csrc) const;

  void string_send(std::string& s, int rdest, int cdest) const;
  void string_recv(std::string& s, int rsrc, int csrc) const;
  void string_bcast(std::string& s, int isrc) const;

  bool operator==(const Context& ctxt) const;

  // MPI communicator for this context. Returns MPI_COMM_NULL if
  // this process is not part of the context
  MPI_Comm comm(void) const;

  // Constructors

  // default global context: construct a single-row global Context
  explicit Context();

  // global Context of size nprow * npcol with column major order
  explicit Context(int nprow, int npcol);

  // construct a Context of size nprow*npcol from the processes
  // in context ctxt starting at process (irstart,icstart)
  explicit Context(const Context &ctxt, int nprow, int npcol,
    int irstart, int icstart);

  ~Context();

  Context(const Context& ctxt);
  Context& operator=(const Context& rhs);

  void print(std::ostream& os) const;
};
std::ostream& operator << ( std::ostream& os, const Context& ctxt );

#endif
