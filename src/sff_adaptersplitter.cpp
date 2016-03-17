/*
 * adapter_dictionary.cpp
 *
 *  Created on: March 17, 2016
 *      Author: Tobias Neumann
 *       Email: tobias.neumann.at@gmail.com
 */

#include "sff_adaptersplitter.h"

/* P R O T O T Y P E S *******************************************************/

// Look up adapter for a single read
void lookup_adapter_per_read(sff_record * record,
		for_ret_dictionary const & dictionary, adapter_map & adapters,
		adapter_record_map & mapping, int const & trim_flag);

void write_common_header(sff_common_header const & h, int const & read_number,
		FILE * out);
void write_record(sff_record * r, sff_common_header const & h, FILE * out);
void write_read_header(sff_read_header * header, FILE * out);
void write_read_data(sff_read_data * data, FILE * out, uint16_t const & nflows,
		uint32_t const & nbases);
void write_padding(FILE * out, int const & header_size);

/* F U N C T I O N S *********************************************************/

// Set up Aho-Corasick tries
for_ret_dictionary set_up_dictionary(adapter_map const & adapters) {

	// Forward dictionary to match a hit at any position
	CSuffixTrie * dictionary = new CSuffixTrie;

	// Reverse dictionary to match only at the starting position
	// Note: Can be easily removed if match at any position should be allowed

	CSuffixTrie * dictionary_rev = new CSuffixTrie;

	for (adapter_map_const_itor i = adapters.begin(); i != adapters.end();
			++i) {
		dictionary->AddString(i->first);
		dictionary_rev->AddString(string(i->first.rbegin(), i->first.rend()));
	}

	dictionary->BuildTreeIndex();
	dictionary_rev->BuildTreeIndex();

	return pair<CSuffixTrie*, CSuffixTrie*>(dictionary, dictionary_rev);

}

// Free allocations
void delete_dictionary(for_ret_dictionary * dictionary) {
	delete dictionary->first;
	dictionary->first = 0;
	delete dictionary->second;
	dictionary->second = 0;
}
