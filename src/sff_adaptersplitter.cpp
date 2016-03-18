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

// Look-up adapters for reads
void map_reads_to_adapters(FILE *sff_fp, adapter_record_map & mapping,
		for_ret_dictionary const & dictionary, adapter_map & adapters,
		uint16_t const & flow_len, int const & numreads, int const & trim_flag) {
	register int cur_pos;

	for (cur_pos = 0; cur_pos < numreads; cur_pos++) {

		sff_record * record = new sff_record;
		sff_read_header * rh = new sff_read_header;
		sff_read_data * rd = new sff_read_data;

		read_sff_read_header(sff_fp, rh);
		read_sff_read_data(sff_fp, rd, flow_len, rh->nbases);

		record->header = rh;
		record->data = rd;

		lookup_adapter_per_read(record, dictionary, adapters, mapping, trim_flag);

	}
}

// Write reads to corresponding adapter file
void split_reads_per_adapter(sff_common_header const & h,
		adapter_record_map & mapping) {

	for (map<string, vector<sff_record *> *>::iterator i = mapping.begin();
			i != mapping.end(); ++i) {

		vector<sff_record *> * current_reads = i->second;
		string filename = i->first + FILE_SUFFIX;

		FILE * out;

		if ((out = fopen(filename.c_str(), "w")) == 0) {
			std::cerr << "[err] Could not open file " << filename
					<< " for reading.\n";
			exit(1);
		}

		write_common_header(h, current_reads->size(), out);

		for (vector<sff_record *>::iterator j = current_reads->begin();
				j != current_reads->end(); ++j) {
			write_record(*j, h, out);
		}

		fclose(out);
	}
}

// Free allocations
void delete_dictionary(for_ret_dictionary * dictionary) {
	delete dictionary->first;
	dictionary->first = 0;
	delete dictionary->second;
	dictionary->second = 0;
}

// Free allocations
void delete_mapping(adapter_record_map * mapping) {
	for (adapter_record_map_itor i = mapping->begin(); i != mapping->end(); ++i) {
		vector<sff_record *> mapped_reads = *i->second;
		for (int cur_read = 0; cur_read < mapped_reads.size(); ++cur_read) {
			free_sff_read_header(mapped_reads[cur_read]->header);
			free_sff_read_data(mapped_reads[cur_read]->data);
			delete mapped_reads[cur_read]->header;
			mapped_reads[cur_read]->header = 0;
			delete mapped_reads[cur_read]->data;
			mapped_reads[cur_read]->data;
			delete mapped_reads[cur_read];
			mapped_reads[cur_read];
		}
		delete i->second;
		i->second = 0;
	}
}

// Look up adapter for a single read
void lookup_adapter_per_read(sff_record * record,
		for_ret_dictionary const & dictionary, adapter_map & adapters,
		adapter_record_map & mapping, int const & trim_flag) {

	// Get bases from clipping
	int left_clip = 0, right_clip = 0, nbases = 0;
	char *bases;

	get_clip_values((*record->header), trim_flag, &left_clip, &right_clip);
	nbases = right_clip - left_clip;
	bases = get_read_bases((*record->data), left_clip, right_clip);

	string seq(bases, nbases);

	CSuffixTrie::DataFound found = dictionary.first->SearchAhoCorasik(seq);

	if (found.sDataFound != "" && found.iFoundPosition == 0) {

		string rev_search(found.sDataFound.rbegin(), found.sDataFound.rend());

		CSuffixTrie::DataFound found_rev = dictionary.second->SearchAhoCorasik(
				rev_search);

		if (found.sDataFound.size() == found_rev.sDataFound.size()) {
			if (mapping[adapters[found.sDataFound]] == 0) {
				mapping[adapters[found.sDataFound]] =
						new vector<sff_record *>();
			}
			mapping[adapters[found.sDataFound]]->push_back(record);
		} else {
			if (mapping[DEFAULT_OUTPUT] == 0) {
				mapping[DEFAULT_OUTPUT] = new vector<sff_record *>();
			}
			mapping[DEFAULT_OUTPUT]->push_back(record);
		}

	} else {
		if (mapping[DEFAULT_OUTPUT] == 0) {
			mapping[DEFAULT_OUTPUT] = new vector<sff_record *>();
		}
		mapping[DEFAULT_OUTPUT]->push_back(record);
	}
}

void write_common_header(sff_common_header const & h, int const & read_number,
		FILE * out) {
	uint32_t bin_read_number = (uint32_t) read_number;

	/* sff files are in big endian notation so adjust appropriately */
	bin_read_number = be32toh(bin_read_number);
	uint32_t be32toh_magic = be32toh(h.magic);
	uint64_t be64toh_index_offset = be64toh(h.index_offset);
	uint32_t be64toh_index_len = be64toh(h.index_len);
	uint16_t be16toh_header_len = be16toh(h.header_len);
	uint16_t be16toh_key_len = be16toh(h.key_len);
	uint16_t be16toh_flow_len = be16toh(h.flow_len);

	fwrite(&be32toh_magic, sizeof(uint32_t), 1, out);
	fwrite(&h.version, sizeof(char), 4, out);
	fwrite(&be64toh_index_offset, sizeof(uint64_t), 1, out);
	fwrite(&be64toh_index_len, sizeof(uint32_t), 1, out);
	fwrite(&bin_read_number, sizeof(uint32_t), 1, out);
	fwrite(&be16toh_header_len, sizeof(uint16_t), 1, out);
	fwrite(&be16toh_key_len, sizeof(uint16_t), 1, out);
	fwrite(&be16toh_flow_len, sizeof(uint16_t), 1, out);
	fwrite(&h.flowgram_format, sizeof(uint8_t), 1, out);
	fwrite(&h.flow, sizeof(char), h.flow_len, out);
	fwrite(&h.key, sizeof(char), h.key_len, out);

	int header_size = sizeof(h.magic) + sizeof(*(h.version)) * 4
			+ sizeof(h.index_offset) + sizeof(h.index_len) + sizeof(h.nreads)
			+ sizeof(h.header_len) + sizeof(h.key_len) + sizeof(h.flow_len)
			+ sizeof(h.flowgram_format) + (sizeof(char) * h.flow_len)
			+ (sizeof(char) * h.key_len);

	if (!(header_size % PADDING_SIZE == 0)) {
		write_padding(out, header_size);
	}
}

void write_record(sff_record * r, sff_common_header const & h, FILE * out) {
	write_read_header(r->header, out);
	write_read_data(r->data, out, h.flow_len, r->header->nbases);
}

void write_read_header(sff_read_header * header, FILE * out) {
	int header_size;

	/* sff files are in big endian notation so adjust appropriately */

	uint16_t be16toh_header_len = be16toh(header->header_len);
	uint16_t be16toh_name_len = be16toh(header->name_len);
	uint32_t be32toh_nbases = be32toh(header->nbases);
	uint16_t be16toh_clip_qual_left = be16toh(header->clip_qual_left);
	uint16_t be16toh_clip_qual_right = be16toh(header->clip_qual_right);
	uint16_t be16toh_clip_adapter_left = be16toh(header->clip_adapter_left);
	uint16_t be16toh_clip_adapter_right = be16toh(header->clip_adapter_right);

	fwrite(&be16toh_header_len, sizeof(uint16_t), 1, out);
	fwrite(&be16toh_name_len, sizeof(uint16_t), 1, out);
	fwrite(&be32toh_nbases, sizeof(uint32_t), 1, out);
	fwrite(&be16toh_clip_qual_left, sizeof(uint16_t), 1, out);
	fwrite(&be16toh_clip_qual_right, sizeof(uint16_t), 1, out);
	fwrite(&be16toh_clip_adapter_left, sizeof(uint16_t), 1, out);
	fwrite(&be16toh_clip_adapter_right, sizeof(uint16_t), 1, out);

	fwrite(header->name, sizeof(char), header->name_len, out);

	/* the section should be a multiple of 8-bytes, if not,
	 it is zero-byte padded to make it so */

	header_size = sizeof(header->header_len) + sizeof(header->name_len)
			+ sizeof(header->nbases) + sizeof(header->clip_qual_left)
			+ sizeof(header->clip_qual_right)
			+ sizeof(header->clip_adapter_left)
			+ sizeof(header->clip_adapter_right)
			+ (sizeof(char) * header->name_len);

	if (!(header_size % PADDING_SIZE == 0)) {
		write_padding(out, header_size);
	}
}

void write_read_data(sff_read_data * data, FILE * out, uint16_t const & nflows,
		uint32_t const & nbases) {

	int data_size;

	/* sff files are in big endian notation so adjust appropriately */

	for (int i = 0; i < nflows; i++) {
		data->flowgram[i] = be16toh(data->flowgram[i]);
	}
	/* read from the sff file and assign to sff_read_data */
	fwrite(data->flowgram, sizeof(uint16_t), (size_t) nflows, out);
	fwrite(data->flow_index, sizeof(uint8_t), (size_t) nbases, out);
	fwrite(data->bases, sizeof(char), (size_t) nbases, out);
	fwrite(data->quality, sizeof(uint8_t), (size_t) nbases, out);

	/* the section should be a multiple of 8-bytes, if not,
	 it is zero-byte padded to make it so */

	data_size = (sizeof(uint16_t) * nflows) // flowgram size
	+ (sizeof(uint8_t) * nbases) // flow_index size
	+ (sizeof(char) * nbases) // bases size
	+ (sizeof(uint8_t) * nbases); // quality size

	if (!(data_size % PADDING_SIZE == 0)) {
		write_padding(out, data_size);
	}
}

void write_padding(FILE * out, int const & header_size) {
	int remainder = PADDING_SIZE - (header_size % PADDING_SIZE);
	uint8_t padding[remainder];
	for (int i = 0; i < remainder; ++i) {
		padding[i] = (uint8_t) 0;
	}
	fwrite(padding, sizeof(uint8_t), remainder, out);
}
