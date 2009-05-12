// $Id$
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

#include <algorithm>
#include "AlignmentElement.h"

using namespace std;

namespace Moses
{
AlignmentElement::AlignmentElement(const ContainerType &alignInfo)
{
	insert_iterator<ContainerType> insertIter( m_collection, m_collection.end() );
	copy(alignInfo.begin(), alignInfo.end(), insertIter);
};

AlignmentElement::AlignmentElement(const vector<AlignmentElementType> &alignInfo)
{
	insert_iterator<ContainerType> insertIter( m_collection, m_collection.end() );
	copy(alignInfo.begin(), alignInfo.end(), insertIter);
};

AlignmentElement::AlignmentElement(const AlignmentElement &alignInfo)
{
	insert_iterator<ContainerType> insertIter( m_collection, m_collection.end() );
	copy(alignInfo.begin(), alignInfo.end(), insertIter);
};

AlignmentElement& AlignmentElement::operator=(const AlignmentElement& alignInfo)
{
	insert_iterator<ContainerType> insertIter( m_collection, m_collection.end() );
	copy(alignInfo.begin(), alignInfo.end(), insertIter);
	
	return *this;
}

void AlignmentElement::Shift(int shift)
{
	ContainerType  newColl;

	ContainerType::const_iterator iter;
	for (iter = m_collection.begin() ; iter != m_collection.end() ; ++iter){
		if (*iter!=-1) newColl.insert(*iter + shift);	
		else newColl.insert(*iter);	
	}
	m_collection = newColl;
}

std::ostream& operator<<(std::ostream& out, const AlignmentElement &alignElement)
{
	const AlignmentElement::ContainerType &elemSet = alignElement.GetCollection();

//	out << "(";
	if (elemSet.size() > 0)
	{
		AlignmentElement::ContainerType::const_iterator iter = elemSet.begin();
		out << *iter;
		for (++iter ; iter != elemSet.end() ; ++iter)
			out << "," << *iter;
	}
//	out << ")";

	return out;
}

void AlignmentElement::SetIntersect(const AlignmentElement &otherElement)
{
	ContainerType newElement;
	set_intersection(m_collection.begin() , m_collection.end()
									,otherElement.begin() , otherElement.end()
									,inserter(newElement , newElement.begin()) );
	m_collection = newElement;
}

void AlignmentElement::SetUniformAlignment(size_t otherPhraseSize)
{
	for (size_t pos = 0 ; pos < otherPhraseSize ; ++pos)
		m_collection.insert(pos);
}

TO_STRING_BODY(AlignmentElement);

}


