/*****************************************************************************

  (c) 2015 - Sun Ruping
  ruping@stanford.edu


g++ changeBamRG.cpp
-I/srv/gsfs0/projects/curtis/ruping/tools/bamtools/include/ -I/srv/gsfs0/projects/curtis/ruping/tools/zlib/current/include/ -I/srv/gsfs0/projects/curtis/ruping/tools/boost/current/include/ 
-L/srv/gsfs0/projects/curtis/ruping/tools/bamtools/lib/ -L/srv/gsfs0/projects/curtis/ruping/tools/zlib/current/lib/ -L/srv/gsfs0/projects/curtis/ruping/tools/boost/current/lib/ 
-lbamtools -lz -Wl,-rpath,/srv/gsfs0/projects/curtis/ruping/tools/bamtools/lib/:/srv/gsfs0/projects/curtis/ruping/tools/boost/current/lib/ -lboost_regex -o changeBamRG
******************************************************************************/

#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamMultiReader.h>

using namespace BamTools;

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <deque>
#include <set>
#include <string>
#include <cstring>
#include <sstream>
#include "changeBamRG.h"
#include <iomanip>
#include "boost/regex.hpp"
using namespace boost;
using namespace std;


inline void splitstring(const string &str, vector<string> &elements, const string &delimiter);
inline string int2str(unsigned int &i);
inline string float2str(float &f);

int main ( int argc, char *argv[] ) { 

  struct parameters *param = 0;
  param = interface(param, argc, argv);


  //bam input and generate index if not yet 
  //-------------------------------------------------------------------------------------------------------+
  // BAM input (file or filenames?)                                                                        |
  //-------------------------------------------------------------------------------------------------------+
  char *fof = param->mapping_f;
  FILE *IN=NULL;
  char linefof[5000];
  int filecount=0;
  vector <string> fnames;

  if (strchr(fof,' ')!=NULL) {
    char *ptr;
    ptr=strtok(fof," ");
    while (ptr!=NULL) {
      fnames.push_back(ptr);
      filecount++;
      ptr=strtok(NULL," ");
    }
  } else {
    IN=fopen(fof,"rt");
    if (IN!=NULL) {
      long linecount=0;
      while (fgets(linefof,5000-1,IN)!=NULL) {
        linecount++;
        if (linefof[0]!='#' && linefof[0]!='\n') {
          char *ptr=strchr(linefof,'\n');
          if (ptr!=NULL && ptr[0]=='\n') {
            ptr[0]='\0';
          }
          FILE *dummy=NULL;
          dummy=fopen(linefof,"rt");
          if (dummy!=NULL) {     // seems to be a file of filenames...
            fclose(dummy);
            fnames.push_back(linefof);
            filecount++;
          } else if (filecount==0 || linecount>=1000-1) {  // seems to be a single file
            fnames.push_back(fof);
            filecount++;
            break;
          }
        }
      }
      fclose(IN);
    }
  }  //file or file name decided and stored in vector "fnames"

  cerr << "the input mapping files are:" << endl;
  vector <string>::iterator fit = fnames.begin();
  for(; fit != fnames.end(); fit++) {
    cerr << *fit << endl;
  }

  //-------------------------------------------------------------------------------------------------------+
  // end of file or filenames                                                                              |
  //-------------------------------------------------------------------------------------------------------+

  // open the BAM file(s)
  BamMultiReader reader;
  reader.Open(fnames);

  // get header & reference information
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();


  // attempt to open BamWriter
  BamWriter writer;
  string outputBam = param->writer;
  if ( outputBam != "" ) {
    if ( !writer.Open(param->writer, header, refs) ) {
      cerr << "Could not open output BAM file" << endl;
      exit(0);
    }
  }


  BamAlignment bam;
  while (reader.GetNextAlignment(bam)) { //change RG 

    bam.EditTag("RG","Z",1)

    writer.SaveAlignment(bam);
    
  }  // read a bam

  return 0;

} //main


inline string int2str(unsigned int &i){
  string s;
  stringstream ss(s);
  ss << i;
  return ss.str();
}


inline string float2str(float &f){
  string s;
  stringstream ss(s);
  ss << f;
  return ss.str();
}


inline void splitstring(const string &str, vector<string> &elements, const string &delimiter) {
  string::size_type lastPos = str.find_first_not_of(delimiter, 0);
  string::size_type pos     = str.find_first_of(delimiter, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    elements.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiter, pos);
    pos = str.find_first_of(delimiter, lastPos);
  }
}
