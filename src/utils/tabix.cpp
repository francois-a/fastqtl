/*
 * This is a C++ wrapper around tabix project which abstracts some of the details of opening and jumping in tabix-indexed files.
 * Author: Erik Garrison erik.garrison@gmail.com
 */

#include "tabix.hpp"

Tabix::Tabix(void) { }

Tabix::Tabix(string& file) {
    filename = file;
    const char* cfilename = file.c_str();
    struct stat stat_tbi,stat_vcf;
    char *fnidx = (char*) calloc(strlen(cfilename) + 5, 1);
    strcat(strcpy(fnidx, cfilename), ".tbi");
    if ( bgzf_is_bgzf(cfilename)!=1 )
    {
        cerr << "[tabix++] was bgzip used to compress this file? " << file << endl;
        free(fnidx);
        exit(1);
    }
    // Common source of errors: new VCF is used with an old index
    stat(fnidx, &stat_tbi);
    stat(cfilename, &stat_vcf);
    if ( stat_vcf.st_mtime > stat_tbi.st_mtime )
    {
        cerr << "[tabix++] the index file is older than the vcf file. Please use '-f' to overwrite or reindex." << endl;
        free(fnidx);
        exit(1);
    }
    free(fnidx);

    if ((t = ti_open(cfilename, 0)) == 0) {
        cerr << "[tabix++] fail to open the data file." << endl;
        exit(1);
    }

    if (ti_lazy_index_load(t) < 0) {
        cerr << "[tabix++] failed to load the index file." << endl;
        exit(1);
    }

    idxconf = ti_get_conf(t->idx);

    // set up the iterator, defaults to the beginning
    iter = ti_query(t, 0, 0, 0);

}

Tabix::~Tabix(void) {
    ti_iter_destroy(iter);
    ti_close(t);
}

void Tabix::getHeader(string& header) {
    header.clear();
    ti_iter_destroy(iter);
    iter = ti_query(t, 0, 0, 0);
    const char* s;
    int len;
    while ((s = ti_read(t, iter, &len)) != 0) {
        if ((int)(*s) != idxconf->meta_char) {
            firstline = string(s); // stash this line
            break;
        } else {
            header += string(s);
            header += "\n";
        }
    }
}

void Tabix::getLastHeader(string & header) {
	header.clear();
	ti_iter_destroy(iter);
	iter = ti_query(t, 0, 0, 0);
	const char* s;
	int len;
	while ((s = ti_read(t, iter, &len)) != 0) {
		if ((int)(*s) != idxconf->meta_char) {
			firstline = string(s); // stash this line
			break;
		} else header = string(s);
	}
}

bool Tabix::setRegion(string region) {
    if (ti_parse_region(t->idx, region.c_str(), &tid, &beg, &end) == 0) {
        firstline.clear();
        ti_iter_destroy(iter);
        iter = ti_queryi(t, tid, beg, end);
        return true;
    } else return false;
}

bool Tabix::getNextLine(string & line) {
    const char* s;
    int len;
    if (!firstline.empty()) {
        line = firstline; // recovers line read if header is parsed
        firstline.clear();
        return true;
    }
    if ((s = ti_read(t, iter, &len)) != 0) {
        line = string(s);
        return true;
    } else return false;
}
