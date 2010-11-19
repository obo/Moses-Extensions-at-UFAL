#include "PhrasePair.h"
#include "Vocabulary.h"

using namespace std;

void PhrasePair::Print( ostream* out, int width ) {
	vector< WORD_ID >::iterator t;
	
	// source
	int sentence_start = m_source_position - m_source_start;
	int source_width = (width-3)/2;
	string source_pre = "";
	string source = "";
	string source_post = "";
	for( int space=0; space<source_width/2; space++ ) source_pre += " ";
	for( char i=0; i<m_source_start; i++ ) {
		source_pre += " " + m_suffixArray->GetWord( sentence_start + i );
	}
	for( char i=m_source_start; i<=m_source_end; i++ ) {
		if (i>m_source_start) source += " ";
		source += m_suffixArray->GetWord( sentence_start + i );
	}
	char source_length = m_suffixArray->GetSentenceLength( m_suffixArray->GetSentence( m_source_position ) );
	for( char i=m_source_end+1; i<source_length; i++ ) {
		if (i>m_source_end+1) source_post += " ";
		source_post += m_suffixArray->GetWord( sentence_start + i );
	}
	for( int space=0; space<source_width/2; space++ ) source_post += " ";
	
	int source_pre_width = (source_width-source.size()-2)/2;
	int source_post_width = (source_width-source.size()-2+1)/2;
	
	if (source.size() > width) {
		source_pre_width = 0;
		source_post_width = 0;
	}
	
	*out << source_pre.substr( source_pre.size()-source_pre_width, source_pre_width ) << " "
	     << source.substr( 0, source_width -2 ) << " "
	     << source_post.substr( 0, source_post_width ) << " | ";
	
	// target
	int target_width = (width-3)/2;
	string target_pre = "";
	string target = "";
	string target_post = "";
	for( int space=0; space<target_width/2; space++ ) target_pre += " ";
	for( char i=0; i<m_target_start; i++ ) {
		target_pre += " " + m_targetCorpus->GetWord( m_sentence_id, i);
	}
	for( char i=m_target_start; i<=m_target_end; i++ ) {
		if (i>m_target_start) target += " ";
		target += m_targetCorpus->GetWord( m_sentence_id, i);
	}
	for( char i=m_target_end+1; i<m_target_length; i++ ) {
		if (i>m_target_end+1) target_post += " ";
		target_post += m_targetCorpus->GetWord( m_sentence_id, i);
	}
	
	int target_pre_width = (target_width-target.size()-2)/2;
	int target_post_width = (target_width-target.size()-2+1)/2;
	
	if (target.size() > width) {
		target_pre_width = 0;
		target_post_width = 0;
	}
	
	*out << target_pre.substr( target_pre.size()-target_pre_width, target_pre_width ) << " "
	     << target.substr( 0, target_width -2 ) << " "
	     << target_post.substr( 0, target_post_width ) << endl;
}

void PhrasePair::PrintTarget( ostream* out ) {
	for( char i=m_target_start; i<=m_target_end; i++ ) {
		if (i>m_target_start) *out << " ";
		*out << m_targetCorpus->GetWord( m_sentence_id, i);
	}
}

void PhrasePair::PrintHTML( ostream* out ) {
	// source
	int sentence_start = m_source_position - m_source_start;
	char source_length = m_suffixArray->GetSentenceLength( m_suffixArray->GetSentence( m_source_position ) );

	*out << "<tr><td align=right class=\"pp_source_left\">";
	for( char i=0; i<m_source_start; i++ ) {
		if (i>0) *out << " ";
		*out << m_suffixArray->GetWord( sentence_start + i );
	}
	*out << "</td><td class=\"pp_source\">";
	for( char i=m_source_start; i<=m_source_end; i++ ) {	
		if (i>m_source_start) *out << " ";
		*out << m_suffixArray->GetWord( sentence_start + i );
	}
	*out << "</td><td class=\"pp_source_right\">";
	for( char i=m_source_end+1; i<source_length; i++ ) {
		if (i>m_source_end+1) *out << " ";
		*out << m_suffixArray->GetWord( sentence_start + i );
	}

	// target
	*out << "</td><td class=\"pp_target_left\">";
	for( char i=0; i<m_target_start; i++ ) {
		if (i>0) *out << " ";
		*out << m_targetCorpus->GetWord( m_sentence_id, i);
	}
	*out << "</td><td class=\"pp_target\">";
	for( char i=m_target_start; i<=m_target_end; i++ ) {
		if (i>m_target_start) *out << " ";
		*out << m_targetCorpus->GetWord( m_sentence_id, i);
	}
	*out << "</td><td class=\"pp_target_right\">";
	for( char i=m_target_end+1; i<m_target_length; i++ ) {
		if (i>m_target_end+1) *out << " ";
		*out << m_targetCorpus->GetWord( m_sentence_id, i);
	}
	*out << "</td></tr>\n";
}

void PhrasePair::PrintClippedHTML( ostream* out, int width ) {
	vector< WORD_ID >::iterator t;
	
	// source
	int sentence_start = m_source_position - m_source_start;
	int source_width = (width+1)/2;
	string source_pre = "";
	string source = "";
	string source_post = "";
	for( char i=0; i<m_source_start; i++ ) {
		source_pre += " " + m_suffixArray->GetWord( sentence_start + i );
	}
	for( char i=m_source_start; i<=m_source_end; i++ ) {
		if (i>m_source_start) source += " ";
		source += m_suffixArray->GetWord( sentence_start + i );
	}
	char source_length = m_suffixArray->GetSentenceLength( m_suffixArray->GetSentence( m_source_position ) );
	for( char i=m_source_end+1; i<source_length; i++ ) {
		if (i>m_source_end+1) source_post += " ";
		source_post += m_suffixArray->GetWord( sentence_start + i );
	}
	int source_pre_width = (source_width-source.size())/2;
	int source_post_width = (source_width-source.size()+1)/2;
	
	if (source.size() > width) {
		source_pre_width = 0;
		source_post_width = 0;
	}
	if (source_pre.size()>source_pre_width)
		source_pre = "..." + source_pre.substr( source_pre.size()-source_pre_width, source_pre_width );
	if (source_post.size() > source_post_width)
		source_post = source_post.substr( 0, source_post_width ) + "...";
	
	*out << "<tr><td class=\"pp_source_left\">"
	     << source_pre 
	     << "</td><td class=\"pp_source\">"
	     << source.substr( 0, source_width -2 )
	     << "</td><td class=\"pp_source_right\">"
	     << source_post
	     << "</td>";
	
	// target
	int target_width = width/2;
	string target_pre = "";
	string target = "";
	string target_post = "";
	for( char i=0; i<m_target_start; i++ ) {
		target_pre += " " + m_targetCorpus->GetWord( m_sentence_id, i);
	}
	for( char i=m_target_start; i<=m_target_end; i++ ) {
		if (i>m_target_start) target += " ";
		target += m_targetCorpus->GetWord( m_sentence_id, i);
	}
	for( char i=m_target_end+1; i<m_target_length; i++ ) {
		if (i>m_target_end+1) target_post += " ";
		target_post += m_targetCorpus->GetWord( m_sentence_id, i);
	}
	
	int target_pre_width = (target_width-target.size())/2;
	int target_post_width = (target_width-target.size()+1)/2;
	
	if (target.size() > width) {
		target_pre_width = 0;
		target_post_width = 0;
	}
	if (target_pre.size() > target_pre_width)
		target_pre = "..." + target_pre.substr( target_pre.size()-target_pre_width, target_pre_width );
	if (target_post.size() > target_post_width)
		target_post = target_post.substr( 0, target_post_width ) + "...";
	
	*out << "<td class=\"pp_target_left\">"
	     << target_pre
	     << "</td><td class=\"pp_target\">"
	     << target.substr( 0, target_width -2 )
	     << "</td><td class=\"pp_target_right\">"
	     << target_post
	     << "</td></tr>"<< endl;
}

