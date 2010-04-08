// $Id: PCNTools.h 3045 2010-04-05 13:07:29Z hieuhoang1972 $

/***********************************************************************
Moses - factored phrase-based language decoder
Copyright (C) 2006 University of Edinburgh

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
***********************************************************************/

#ifndef moses_PCNTools
#define moses_PCNTools

#include <vector>
#include <string>
#include <utility>
#include <cstdlib>

/** A couple of utilities to read .pcn files. A python-compatible format
  * for encoding confusion networks.
  */
namespace PCN {

  typedef std::pair<std::pair<std::string, std::vector<float> >, size_t> CNAlt;
  typedef std::vector<CNAlt> CNCol;
  typedef std::vector<CNCol> CN;

  /** Given a string ((('foo',0.1),('bar',0.9)),...) representation of a
    * confusion net in PCN format, return a CN object 
    */
  CN parsePCN(const std::string& in);
  
};

#endif
