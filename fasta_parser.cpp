#include <iostream>
#include <getopt.h>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <assert.h>
#include <stdio.h>
//#include <regex>

using namespace std;

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1"
#endif


char nst_dent4_table[4]={'A', 'C', 'G', 'T'};
uint64_t pac_step = (64 * 1024 * 1024);

typedef struct condOpt{
    string input;
    string output;
}CondOpt;

typedef struct ann_t{
    string contig;
    long start;
    long length;
    map<uint64_t,uint64_t> N_region;
}Ann_t;

typedef struct amb_t{;
    long start;
    long length;
    char C;
}Amb_t;

CondOpt opt;
//ann amb file info vec
vector<ann_t > ann_info;
vector<amb_t > amb_info;

void usage(){

    fprintf(stderr, "\n");
    fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
    fprintf(stderr, "Build at: %s  %s, GCC version: %d.%d.%d \n", __DATE__, __TIME__, __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
    fprintf(stderr, "Usage:   fasta_parser [-i input.fa]\n\n");
    fprintf(stderr, "options: \n");
    fprintf(stderr, "         i,input          pac, ann, amb input file. \n");
    fprintf(stderr, "         o,output         fastA output file.\n");
    fprintf(stderr, "         h,help           help\n");
    fprintf(stderr, "\n");

}

//split string with 'c' and return vector size.
size_t split_string(const std::string &s, std::vector<std::string> &v, const std::string &c)
{
    std::string::size_type pos1, pos2;
    pos2 = s.find(c);
    pos1 = 0;
    while(std::string::npos != pos2)
    {
        v.push_back(s.substr(pos1, pos2-pos1));

        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if(pos1 != s.length())
        v.push_back(s.substr(pos1));

    return v.size();
}

int main(int argc, char* argv[]){

    // get opt
    const char *shortOptions = "hi:o:";
    const struct option longOptions[] =
            {
                    { "help", 1, NULL, 'h' },
                    { "input", 1, NULL, 'i' },
                    { "output", 1, NULL, 'o' },
            };

    if(argc<2) usage(),exit(0);
    int nextOpt;
    while(-1 != (nextOpt = getopt_long(argc, argv, shortOptions, longOptions, NULL))) {
        switch (nextOpt) {
            case 'i':
                opt.input = string(optarg);
                opt.output = string(optarg)+".dump";
                break;
            case 'o':
                opt.output = string(optarg);
                break;
            case 'h':
                usage();
                break;
        }
    }

    // open file.
    ifstream input_ann(opt.input+".ann");
    if(!input_ann.is_open()){
        cerr<<"[Error]: Fail open file:"<<opt.input+".ann"<<endl;
        exit(1);
    }
    ifstream input_amb(opt.input+".amb");
    if(!input_amb.is_open()){
        cerr<<"[Error]: Fail open file:"<<opt.input+".amb"<<endl;
        exit(1);
    }
    ifstream input_pac(opt.input+".pac", ios::in|ios::binary);
    if(!input_pac.is_open()){
        cerr<<"[Error]: Fail open file:"<<opt.input+".pac"<<endl;
        exit(1);
    }
    ofstream output_fatsa(opt.output, ios::out|ios::trunc);
    if(!output_fatsa.is_open()){
        cerr<<"[Error]: Fail open file:"<<opt.output<<endl;
        exit(1);
    }


    char input_buff[512];
    //parse ann info
    input_ann.getline(input_buff, 512);
    while(!input_ann.eof()){
        vector<string > ann_line;
        input_ann.getline(input_buff, 512);
        if(split_string(string(input_buff), ann_line, " ")){   //scan first line
            input_ann.getline(input_buff, 512);
            if(split_string(string(input_buff), ann_line, " ")) {   //scan second line
                Ann_t ann = {
                        .contig = ann_line[1],
                        .start = stol(ann_line[3]),
                        .length = stol(ann_line[4]),
                };
                ann_info.emplace_back(ann);
            }
        }
    }

    //parse amb info
    input_amb.getline(input_buff, 512);
    while(!input_amb.eof()){
        vector<string > amb_line;
        input_amb.getline(input_buff, 512);
        if(split_string(string(input_buff), amb_line, " ")){   //scan first line
            Amb_t amb = {
                    .start = stol(amb_line[0]),
                    .length = stol(amb_line[1]),
                    .C = amb_line[2].c_str()[0],
            };
            amb_info.emplace_back(amb);
            //insert amb info into ann info
            for(ann_t &ann : ann_info){
                if(amb.start>=ann.start && amb.start<ann.start+ann.length) {
                    ann.N_region.insert(make_pair(amb.start-ann.start, amb.length));
                    break;
                }
            }
        }
    }

    //print parser messages
    cout<<"[Message]: success parse ann file, size:"<< ann_info.size() <<endl;

//    for(ann_t ann : ann_info){
//        cout<< ann.contig<<":"<<ann.start<<"~"<<ann.start+ann.length<<endl;
//        for(auto n: ann.N_region){
//            cout<<"\t"<<n.first<<"+"<<n.second<<endl;
//        }
//    }

    //read pac file
    uint64_t id=1;
    for(ann_t ann : ann_info){
        cout<<"[Message]: processe ann "<< id++ <<"/"<< ann_info.size() << "\tname: "<< ann.contig << "\tlength: "<< ann.length <<endl;
        output_fatsa << '>' + ann.contig << endl;

        if(ann.length==0) continue; //skip this ann
        uint64_t pac_start_x = ann.start/4;
        uint64_t pac_start_y = ann.start%4; //pac offset
        //if pac_start_y!=0, pac need read one more byte
        uint64_t pac_length = pac_start_y?(ann.length-1)/4+2:(ann.length-1)/4+1;

        char* pac_buff = (char*) calloc(pac_length, sizeof(char));
        char* dump_buff = (char*) calloc(1+pac_length*4, sizeof(char));
        assert(pac_buff!=nullptr);
        assert(dump_buff!= nullptr);
        input_pac.seekg(pac_start_x);
        input_pac.read(pac_buff, pac_length );


        for(uint64_t i=0; i<pac_length; i++){
            dump_buff[4*i+0] = nst_dent4_table[(pac_buff[i]>>6)&0b11] ;
            dump_buff[4*i+1] = nst_dent4_table[(pac_buff[i]>>4)&0b11] ;
            dump_buff[4*i+2] = nst_dent4_table[(pac_buff[i]>>2)&0b11] ;
            dump_buff[4*i+3] = nst_dent4_table[(pac_buff[i]>>0)&0b11] ;
        }
        dump_buff[ann.length+pac_start_y]=0;    //cut the input_buff tail
        //print N region
        for(auto n: ann.N_region){
            assert(pac_start_y+n.first+n.second<=pac_length*4);
            memset(dump_buff+pac_start_y+n.first,'N', n.second);
        }

        if(1){
            //output fa with strip 50
            char out_buff[51]={0};
            uint64_t i=0;
            for(i=0; i<ann.length-50; i+=50){
                memcpy(out_buff, dump_buff+pac_start_y +i, 50);
                output_fatsa<<out_buff<<endl;
            }
            uint8_t lazy_dump_buff_size = ann.length+pac_start_y-i+1;
            memcpy(out_buff, dump_buff+pac_start_y +i, lazy_dump_buff_size);
            output_fatsa<<out_buff<<endl;
        } else {
            //output fa without strip
            output_fatsa<<dump_buff+pac_start_y<<endl;
            output_fatsa.write(dump_buff+pac_start_y, ann.length);
            output_fatsa<<endl;
        }
        free(pac_buff);
        free(dump_buff);
    }

    output_fatsa.close();
    input_ann.close();
    input_amb.close();
    input_pac.close();

}