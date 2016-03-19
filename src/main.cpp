//============================================================================
// Name        : sff-adaptersplitter.cpp
// Author      : Tobias Neumann
// Version     : 1.0
// Email       : tobias.neumann.at@gmail.com
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#include "sff_adaptersplitter.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::stringstream;

/* D E F I N E S *************************************************************/
#define PRG_NAME 	"sff-adaptersplitter"
#define VERSION 	"1.0"
#define ARG_COUNT 	3

/* P R O T O T Y P E S *******************************************************/

void process_options(int argc, char *argv[]);
void processing_message();
void shutdown_message(uint32_t const & read_number);
void help_message();
adapter_map parse_adapters();

/* G L O B A L S *************************************************************/

// Trimming set on per default
int trim_flag = 1;

sff_common_header h;
string adapter_filename = "";
string sff_filename = "";

/* F U N C T I O N S *********************************************************/

int main(int argc, char *argv[]) {

	process_options(argc, argv);
	processing_message();
	parse_adapters();

	FILE * sff_fp;
	adapter_map adapters;
	for_ret_dictionary dictionary;
	adapter_record_map mapping;

	// Parse adapter file and set up Aho-Corasick tries

	adapters = parse_adapters();
	dictionary = set_up_dictionary(adapters);

	// Open sff-file

	if ((sff_fp = std::fopen(sff_filename.c_str(), "r")) == NULL) {
		cerr << "[err] Could not open file " << sff_filename
				<< " for reading.\n";
		exit(1);
	}

	// Common header
	read_sff_common_header(sff_fp, &h);
	verify_sff_common_header(PRG_NAME, VERSION, &h);

	int const numreads = (int) h.nreads;

	// Lookup adapters for reads
	map_reads_to_adapters(sff_fp, mapping, dictionary, adapters, h.flow_len,
			numreads, trim_flag);

	// Write reads to corresponding adapter file
	split_reads_per_adapter(h, mapping);

	// Free allocations
	delete_dictionary(&dictionary);
	delete_mapping(&mapping);

	shutdown_message(h.nreads);

	free_sff_common_header(&h);

	fclose(sff_fp);

	return 0;
}

void process_options(int argc, char *argv[]) {
	if (argc == 1) {
		help_message();
		exit(0);
	} else if (argc == ARG_COUNT) {
		adapter_filename = argv[1];
		sff_filename = argv[2];
	} else if (argc == ARG_COUNT + 1) {
		adapter_filename = argv[1];
		sff_filename = argv[2];
		trim_flag = atoi(argv[3]);
	} else {
		cerr << "Error: Illegal parameters! (see usage)\n";
		exit(1);

	}
}

void processing_message() {
	string trimming;
	trim_flag ? trimming = "on" : trimming = "off";

	cout << "[Message]\tStarting up sff_splitByAdapter:" << endl
			<< "[Message]\tRunning on:\t" << sff_filename << endl
			<< "[Message]\tAdapters from:\t" << adapter_filename << endl
			<< "[Message]\tTrimming before matching adapters:\t" << trimming
			<< endl;
}

void shutdown_message(uint32_t const & read_number) {
	cout << "[Message]\tProcessing finished (" << read_number
			<< " reads were processed)." << endl;
}

void help_message() {
	cout
			<< "Usage: "
			<< PRG_NAME << " <adapter file> <sff file> | <trimming 0|1 (default: 1)>\n";
}

adapter_map parse_adapters() {

	adapter_map adapters;
	ifstream adapter_stream(adapter_filename.c_str());

	if (adapter_stream.is_open()) {
		while (adapter_stream.good()) {
			string line;
			getline(adapter_stream, line);

			if (line.size() != 0) {
				stringstream linestream(line);
				string adapter, seq;

				linestream >> adapter >> seq;

				if (adapter.size() == 0 || seq.size() == 0) {
					cerr
							<< "[Error] Error in adapter file processing! (check format)"
							<< endl;
					exit(1);
				}
				adapters[seq] = adapter;
			}
		}
		adapter_stream.close();
	} else {
		cerr << "[Error] Could not open adapter file " << adapter_filename
				<< "!\n";
		exit(1);
	}

	return adapters;
}
