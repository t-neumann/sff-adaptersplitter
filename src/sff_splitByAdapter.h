/*
 * adapter_dictionary.h
 *
 *  Created on: March 15, 2016
 *      Author: Tobias Neumann
 *       Email: tobias.neumann.at@gmail.com
 */

#ifndef ADAPTER_DICTIONARY_H_
#define ADAPTER_DICTIONARY_H_

#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "./sff/sff.h"
#include "./AhoCorasick/StdAfx.h"
#include "./AhoCorasick/SuffixTrie.h"

using std::map;
using std::string;
using std::vector;
using std::pair;
using std::cerr;

// ---------------- DEFINES  ---------------------------------------

#define DEFAULT_OUTPUT "no_match"
#define FILE_SUFFIX ".sff"

// ----------------- TYPES  ----------------------------------------

typedef map<string, string> adapter_map;
typedef map<string, string>::iterator adapter_map_itor;
typedef map<string, string>::const_iterator adapter_map_const_itor;

typedef pair<CSuffixTrie *, CSuffixTrie*> for_ret_dictionary;

typedef struct {
	sff_read_header * header;
	sff_read_data * data;
} sff_record;

typedef map<string, vector<sff_record *> *> adapter_record_map;
typedef map<string, vector<sff_record *> *>::iterator adapter_record_map_itor;

// ----------------- INTERFACES  -----------------------------------

// Set up Aho-Corasick tries
for_ret_dictionary set_up_dictionary(adapter_map const & adapters);
// Look-up adapters for reads
void map_reads_to_adapters(FILE *sff_fp, adapter_record_map & mapping,
		for_ret_dictionary const & dictionary, adapter_map & adapters,
		uint16_t const & flow_len, int const & numreads, int const & trim_flag);
// Write reads to corresponding adapter file
void split_reads_per_adapter(sff_common_header const & h,
		adapter_record_map & mapping);

// Free allocations
void delete_dictionary(for_ret_dictionary * dictionary);
void delete_mapping(adapter_record_map * mapping);

#endif /* ADAPTER_DICTIONARY_H_ */
