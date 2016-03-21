/*
 * This is a C++ wrapper around tabix project which abstracts some of the details of opening and jumping in tabix-indexed files.
 * Author: Erik Garrison erik.garrison@gmail.com
 */

#ifndef _TABIX_HPP
#define _TABIX_HPP

#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include <bgzf.h>
#include <tabix.h>
#include <iostream>


using namespace std;

class Tabix {

    tabix_t *t;
    ti_iter_t iter;
    const ti_conf_t *idxconf;
    int tid, beg, end;
    string firstline;

public:

    string filename;

    Tabix(void);
    Tabix(string& file);
    ~Tabix(void);

    void getHeader(string& header);
    void getLastHeader(string & header);
    bool setRegion(string region);
    bool getNextLine(string& line);
};

#endif
