///////////////////////////////////////////////////////////////////////////////
// maxprotein.hh
//
// Compute the set of foods that maximizes protein, within a calorie budget,
// with the greedy method or exhaustive search.
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

// Simple structure for a single protein
struct Protein {
    Protein() {
        description = "";
        sequence = "";
    }
    Protein(std::string desc, std::string seq) {
        description = desc;
        sequence = seq;
    }
    std::string		description;
    std::string 	sequence;
};

// Alias for a vector of shared pointers to Protein objects.
typedef std::vector<std::shared_ptr<Protein>> ProteinVector;


// -------------------------------------------------------------------------
// Load all the proteins from a standard FASTA format file with one line
// per sequence (multi-line sequences are not allowed).
// Returns false on I/O error.
// -------------------------------------------------------------------------
bool load_proteins(ProteinVector & proteins, const std::string& path)
{
    //std::cout << "Loading proteins from [" << path << "]" << std::endl;
    proteins.clear();
    std::ifstream ifs(path.c_str());
    if (!ifs.is_open() || !ifs.good()) {
        std::cout << "Failed to open [" << path << "]" << std::endl;
        return false;
    }
    int proteinsLoaded = 0;
    bool have_description = false;
    std::shared_ptr<Protein> newProtein = nullptr;
    while (!ifs.eof()) {
        std::string lineBuffer;
        std::getline(ifs, lineBuffer);
        if (ifs.eof()) {
            break;
        }
        if (lineBuffer.size() == 0) {
            continue;
        }
        if (lineBuffer[0] == '>') {
            newProtein = std::shared_ptr<Protein>(new Protein);
            newProtein->description = lineBuffer.substr(1);
            have_description = true;
        } else if (have_description) {
            newProtein->sequence = lineBuffer;
            proteins.push_back(newProtein);
            proteinsLoaded++;
            have_description = false;
        }
    }
    
    ifs.close();
    //std::cout << "Loaded " << proteinsLoaded << " proteins from [" << path << "]" << std::endl;
    
    return true;
}


// -------------------------------------------------------------------------
int dynamicprogramming_longest_common_subsequence(const std::string & string1,
                                                  const std::string & string2)
{
    int n = string1.size(), m = string2.size();
    if (string1.empty() || string2.empty())
        return 0;
    
    std::vector< std::vector<int> > tbl(n + 1, std::vector<int>(m + 1, 0));
   
    int i = 0, j = 0;

    std::vector< std::vector<int> > up(i, j);
    std::vector< std::vector<int> > left(i, j);
    std::vector< std::vector<int> > diag(i, j);

    for(i = 1; i <= n; i++) {
        for(j = 1; j <= m; j++) {
            up[i][j] = tbl[i - 1][j];
	    left[i][j] = tbl[i][j - 1];
	    diag = tbl[i - 1][j - 1];

	    if(string1[i - 1] == string2[j - 1])
	    	diag[i][j] += 1;	
	
	    tbl[i][j] = fmax(up[i][j], left[i][j]);
	    tbl[i][j] = fmax(left[i][j], diag[i][j]);
	    tbl[i][j] = fmax(up[i][j], diag[i][j]);
        }
    }
    
    return tbl[n][m];
}

// -------------------------------------------------------------------------
std::unique_ptr<std::vector<std::string>> generate_all_subsequences(const std::string & sequence)
{
    std::unique_ptr<std::vector<std::string>> subsequences(nullptr);
    subsequences = std::unique_ptr<std::vector<std::string>> (new std::vector<std::string>);

    unsigned int pow_set_size = pow(2, sequence.size());
    short counter, j;
    std::string subsequence = "";

    for(counter = 0; counter < pow_set_size;  counter++) {
        for(j = 0; j <  sequence.size(); j++) {
            short c = (counter>>j)&1 ;

            if(c == 1)
                subsequence += sequence[j];
        }
        if(find(subsequences -> begin(), subsequences -> end(), subsequence) == subsequences -> end())
            subsequences -> push_back(subsequence);
    }
    
    return subsequences;
   }

// -------------------------------------------------------------------------
int exhaustive_longest_common_subsequence(const std::string & string1, const std::string & string2)

{
    std::unique_ptr<std::vector<std::string>>all_subseqs1(nullptr);
    all_subseqs1= std::unique_ptr<std::vector<std::string>> (new std::vector<std::string>);
    all_subseqs1 = generate_all_subsequences(string1);

    std::unique_ptr<std::vector<std::string>>all_subseqs2(nullptr);
    all_subseqs2= std::unique_ptr<std::vector<std::string>> (new std::vector<std::string>);
    all_subseqs2 = generate_all_subsequences(string2);

    int best_score = 0, s1, s2;

    for(s1 = 0; s1 < all_subseqs1->size(); ++s1) {
        for(s2 = 0; s2 < all_subseqs2->size(); ++s2) {
	    if(string1 == string2 && string1.length() > best_score) 
                best_score = string1.length();
        }
    }

    return best_score;
}

// -------------------------------------------------------------------------
std::shared_ptr<Protein> exhaustive_best_match(ProteinVector & proteins, const std::string & string1)
{
    int score = 0;
    int best_score = 0;
    int best_i = 0;

    std::shared_ptr<Protein> best_protein = nullptr;
    best_protein = std::shared_ptr<Protein> (new Protein);

    for(int i = 0; i < proteins.size(); i++) {    
        score = exhaustive_longest_common_subsequence(proteins.at(i).get() -> sequence,string1);
        if(score > best_score) {
            best_score = score;
            best_i = i;
        }
    }

    best_protein = proteins[best_i];
    
    return best_protein;
}

// -------------------------------------------------------------------------
std::shared_ptr<Protein> dynamicprogramming_best_match(ProteinVector & proteins, const std::string & string1)
{ 
    int best_i = 0, best_score = 0, score;

    std::shared_ptr<Protein> best_protein = nullptr;
    best_protein =  std::shared_ptr<Protein> (new Protein);
    
    for(int i = 0; i < proteins.size(); i++) {
        score = dynamicprogramming_longest_common_subsequence(proteins.at(i).get() -> sequence, string1);
        if(score > best_score) {
            best_score = score;
            best_i = i;
        }
    }
    
    best_protein = proteins[best_i];

    return best_protein;
    
}