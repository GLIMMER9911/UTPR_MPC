// Copyright 2018 The MathWorks, Inc.

#include <vector>
#include <cstring>
#include <string>
#include <array>
#include "statcal_util.hpp"


// [int len][double][double]...[double]
// [int -len][char][char]...[char]

int decode_header(const char *str)
{
    int len;
    std::memcpy(&len, str, sizeof(int));
    return len;
}

void encode_double_data(const std::array<double, MAX_SIZE> & data, const size_t size, std::string & ec)
{
    size_t rlen = sizeof(double)*size;

    ec.resize(sizeof(int)+rlen, 0);
    std::memcpy(&ec[0], &size, sizeof(int));
    std::memcpy(&ec[sizeof(int)], &data[0], rlen);
}

void encode_str_data(const char *data, std::string & ec)
{
    int len = strlen(data);
    size_t rlen = len*sizeof(char);

    ec.resize(rlen+sizeof(int), 0);

    len = -len;
    std::memcpy(&ec[0], &len, sizeof(int));
    std::memcpy(&ec[sizeof(int)], data, rlen);
}


void decode_double_data(const int len, const char *str, std::vector<double> & data)
{
    data.resize(len, 0);
    std::memcpy(&data[0], str+sizeof(int), len*sizeof(double)); 
}

void decode_str_data(const char *str, std::string &s)
{
    s = str+sizeof(int);
}
