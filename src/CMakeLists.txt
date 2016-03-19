cmake_minimum_required(VERSION 2.8)
project(SffAdapterSplitter)

ADD_DEFINITIONS(-g -O2)

add_executable(sff-adaptersplitter 	main.cpp
				sff_adaptersplitter.cpp
				AhoCorasick/StdAfx.cpp
				AhoCorasick/SuffixTrie.cpp
				sff/sff.c
	      )


# Installation of the program
INSTALL(TARGETS sff-adaptersplitter 
        DESTINATION "bin"
        PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
