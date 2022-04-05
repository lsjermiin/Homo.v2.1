/////////////////////////////////////////////////////////////////////////////////
// Program name       : homo.cpp
//
// Version            : See VERSION under DECLARATION OF EXTERNAL VARIABLES
//
// Author             : Lars S Jermiin
//
// Institutions       : Australian National University
//                      Research School of Biology
//                      Acton, ACT 2601, Australia
//
//                      Univerity College Dublin
//                      School of Biology & Environmental Science
//                      Belfield, Dublin 4, Ireland
//
// Emails             : lars.jermiin [at] anu.edu.au
//                      lars.jermiin [at] ucd.ie
//
// URL                : https://github.com/lsjermiin/Homo2.1
//
// Date begun         : 14 April, 2019
//
// Date modified      : 25 March, 2022
//
// Copyright          : Copyright © 2019-22 Lars Sommer Jermiin.
//                      All rights reserved.
//
// Responsibility     : The copyright holder takes no legal responsibility for
//                      the correctness of results obtained using this program.
//
// Summary            : Homo is designed to conduct the matched-pairs test of
//                      symmetry (Bowker 1948) for alignments of:
//
//                        1   Nucleotides (A|C|G|T) (4 states)
//                        2   Nucleotides recoded CTR = (C|T|AG) (3 states)
//                        3   Nucleotides recoded AGY = (A|G|CT) (3 states)
//                        4   Nucleotides recoded ATS = (A|T|CG) (3 states)
//                        5   Nucleotides recoded CGW = (C|G|AT) (3 states)
//                        6   Nucleotides recoded ACK = (A|C|GT) (3 states)
//                        7   Nucleotides recoded GTM = (G|T|AC) (3 states)
//                        8   Nucleotides recoded KM = (GT|AC)   (2 states)
//                        9   Nucleotides recoded RY = (AG|CT)   (2 states)
//                       10   Nucleotides recoded SW = (GC|AT)   (2 states)
//                       11   Nucleotides recoded AB = (A|CGT)   (2 states)
//                       12   Nucleotides recoded CD = (C|AGT)   (2 states)
//                       13   Nucleotides recoded GH = (G|ACT)   (2 states)
//                       14   Nucleotides recoded TV = (T|ACG)   (2 states)
//                       15   Di-nucleotides (AA|AC|...|TG|TT)  (16 states)
//                       16   Di-nucleotides 1st position (A|C|G|T) (4 states)
//                       17   Di-nucleotides 2nd position (A|C|G|T) (4 states)
//                       18   Codons (AAA|AAC|...|TTG|TTT) (64 states)
//                       19   Codons 1st + 2nd positions (AA|AC|...|TG|TT) (16 states)
//                       20   Codons 1st + 3rd positions (AA|AC|...|TG|TT) (16 states)
//                       21   Codons 2nd + 3rd positions (AA|AC|...|TG|TT) (16 states)
//                       22   Codons 1st + 2nd positions (A|C|G|T) (4 states)
//                       23   Codons 1st + 3rd positions (A|C|G|T) (4 states)
//                       24   Codons 2nd + 3rd positions (A|C|G|T) (4 states)
//                       25   Codons 1st position (A|C|G|T) (4 states)
//                       26   Codons 2nd position (A|C|G|T) (4 states)
//                       27   Codons 3rd position (A|C|G|T) (4 states)
//                       28   10-state genotype data (A|C|G|T|K|M|R|Y|S|W) (10 states)
//                       29   14-state genotype data (A|C|G|T|K|M|R|Y|S|W|B|D|H|V) (14 states)
//                       30   Amino acids (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C) (20 states)
//                       31   Recoded amino acids Dayhoff-6 = (AGPST|DENQ|HKR|MIVL|WFY|C) (6 states)
//
//                      For each pair of sequences, Homo generates a P value (i.e., the
//                      probability of getting a chi-square distributed random variable
//                      that equals or exceeds the observed test statistic.
//
//                      The multiple comparison problem is dealt with by controlling the
//                      the family-wise error rate and the false discovery rate.
//
//                      For each type of data, Homo also computes four distances:
//                        *   p distance
//                        *   Compositional distance (using Bowker's test statistic)
//                        *   Compositional distance (using Euclidean's metric - full symmetry)
//                        *   Compositional distance (using Euclidean's metric - marginal symmetry)
//
//                      The input for each distance is the data in a divergence
//                      matrix.
//
// Data               : Sequences must be stored in the FASTA format.
//
// Processing         : Characters are converted to integers to speed up the
//                      program.
//
// Nucleotides        : Alphabet: [A,C.G,T/U,-] = [0,1,2,3,4].
//
//                      Ambiguous characters (i.e., ?, N, B, D, H, K, M, R, S,
//                      V, W and Y) are treated as if they were alignment gaps
//                      (-) (i.e., as missing data).
//
//                      10-state genotypes : Alphabet: [A,C,G,K,M,R,S,T/U,W,Y,-] =
//                      [0,1,2,3,4,5,6,7,8,9,10].
//
//                      Ambiguous characters (i.e., ?, N, B, D, G and V) are
//                      treated as if they were alignment gaps (-) (i.e., as
//                      missing data).
//
//                      14-state genotypes : Alphabet: [A,C,G,T/U,K,M,R,S,W,Y,B,D,H,V,-] =
//                      [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14].
//
//                      Ambiguous characters (i.e., ? and N) are treated as if
//                      they were alignment gaps (-) (i.e., as missing data).
//
// Amino acids        : Alphabet: [A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,-] =
//                      [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
//
//                      Ambiguous characters (i.e., ?, X and Z) are treated as
//                      if they were alignment gaps (-) (i.e., as missing data).
//
// Manuscript         : Jermiin et al. (2022).
//                      On making sense of multiple P values from matched-pairs
//                      tests of homogeneity for pairs of homologous sequences.
//                      Syst. Biol. (in prep.)
//
// References         : Benjamini Y., Yekutieli D. (2001). The control of false
//                      discovery rate in multiple testing under dependency. Ann.
//                      Statist. 29, 1165-1188
//
//                      Bonferroni C. E. (1936). Teoria statistica delle classi e
//                      calcolo delle probabilità. Pubblicazioni del R Istituto
//                      Superiore di Scienze Economiche e Commerciali di Firenze
//                      8, 3-62.
//
//                      Bowker A. H. (1948). A test of symmetry in contingency
//                      tables. J. Am. Stat. Assoc. 43, 572-574.
//
//                      Holm S. (1979). A simple sequentially rejective multiple
//                      test procedure. Scand. J. Stat. 6, 65-70.
//
// Revisions/debugs   : 19 Oct 2021
//                      Scaled d_cfs (i.e. d_cfs = d_cfs/sum_dm) so the value is
//                      given per site.
//
//                      5 Mar 2022
//                      Corrected compositional distance d_Bowker
//
//                      27 Feb 2022 - 25 Mar 2022
//                      Major redesign of the main function to include methods to
//                      control the family-wise error rate (FWER) as well as the
//                      false discovery rate (FDR)
/////////////////////////////////////////////////////////////////////////////////

#include <cctype>
#include <cmath>
#include <limits>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#define SQR(a) ((a) * (a))

using namespace std;
using std::string;
using std::vector;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

// DECLARATION OF EXTERNAL VARIABLES
// The following variables are declared here as they are needed in different functions

const unsigned TWO(2);          // for 2-state alphabet (recoded DNA)
const unsigned THREE(3);        // for 3-state alphabet (recoded DNA)
const unsigned FOUR(4);         // for 4-state alphabet (DNA)
const unsigned SIX(6);          // for 6-state alphabet (recoded amino acids)
const unsigned TEN(10);         // for 10-state alphabet (genotype data)
const unsigned FOURTEEN(14);    // for 14-state alphabet (genotype data)
const unsigned SIXTEEN(16);     // for 16-state alphabet (subset of codons)
const unsigned TWENTY(20);      // for 20-state alphabet (amino acids)
const unsigned SIXTYFOUR(64);   // for 64-state alphabet (codons)
const unsigned max_array(65);   // corresponding to codons & gaps
const double VERSION(2.1);      // version number
vector<string> taxon;           // 2D container for sequence names
vector<vector<int> > alignment; // 2D container for sequence data


// DECLARATION OF FUNCTIONS
// This function translates a string of characters into a vector of integers
vector<int> Translator(unsigned datatype, string seq) {
    int unit; // integer for singlet, duplet or triplet (codon)
    string duplet(""), triplet(""); // strings for dinucleotides and codons
    vector<int> seq_data;
    
    switch (datatype) {
        case 1: // Nucleotides (A|C|G|T)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 2: // Nucleotides (C|T|R)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'C': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'A': seq_data.push_back(2); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'R': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 3: // Nucleotides (A|G|Y)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(2); break;
                    case 'U': seq_data.push_back(2); break;
                    case 'Y': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 4: // Nucleotides (A|T|S)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'C': seq_data.push_back(2); break;
                    case 'S': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 5: // Nucleotides (C|G|W)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'C': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'A': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(2); break;
                    case 'U': seq_data.push_back(2); break;
                    case 'W': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 6: // Nucleotides (A|C|K)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(2); break;
                    case 'U': seq_data.push_back(2); break;
                    case 'K': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 7: // Nucleotides (G|T|M)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'G': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'A': seq_data.push_back(2); break;
                    case 'C': seq_data.push_back(2); break;
                    case 'M': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 8: // Nucleotides (K|M)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'G': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(0); break;
                    case 'U': seq_data.push_back(0); break;
                    case 'K': seq_data.push_back(0); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'M': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 9: // Nucleotides (R|Y)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(0); break;
                    case 'R': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'Y': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 10: // Nucleotides (S|W)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'C': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(0); break;
                    case 'R': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'Y': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 11: // Nucleotides (A|B)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'B': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 12: // Nucleotides (C|D)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'C': seq_data.push_back(0); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'D': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 13: // Nucleotides (G|H)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'G': seq_data.push_back(0); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'H': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 14: // Nucleotides (T|V)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'T': seq_data.push_back(0); break;
                    case 'U': seq_data.push_back(0); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'V': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 15: // Di-nucleotides (AA|AC|..|TG|TT)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 2) {
                for (string::size_type j = i; j != i + 2; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': duplet.push_back('0'); break;
                        case 'C': duplet.push_back('1'); break;
                        case 'G': duplet.push_back('2'); break;
                        case 'T': duplet.push_back('3'); break;
                        case 'U': duplet.push_back('3'); break;
                        default : duplet.push_back('4'); break;
                    }
                }
                unit = stoi(duplet);
                duplet.clear();
                switch (unit) {
                    case 00: seq_data.push_back(0);  break; // AA
                    case 01: seq_data.push_back(1);  break; // AC
                    case 02: seq_data.push_back(2);  break; // AG
                    case 03: seq_data.push_back(3);  break; // AT
                    case 10: seq_data.push_back(4);  break; // CA
                    case 11: seq_data.push_back(5);  break; // CC
                    case 12: seq_data.push_back(6);  break; // CG
                    case 13: seq_data.push_back(7);  break; // CT
                    case 20: seq_data.push_back(8);  break; // GA
                    case 21: seq_data.push_back(9);  break; // GC
                    case 22: seq_data.push_back(10); break; // GG
                    case 23: seq_data.push_back(11); break; // GT
                    case 30: seq_data.push_back(12); break; // TA
                    case 31: seq_data.push_back(13); break; // TC
                    case 32: seq_data.push_back(14); break; // TG
                    case 33: seq_data.push_back(15); break; // TT
                    default: seq_data.push_back(16); break; // In case of other characters
                }
            }
            break;
        case 16: // Di-nucleotides 1st position (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 2) {
                switch (toupper(seq[i])) {
                    case 'A' : seq_data.push_back(0); break; // 1st A
                    case 'C' : seq_data.push_back(1); break; // 1st C
                    case 'G' : seq_data.push_back(2); break; // 1st G
                    case 'T' : seq_data.push_back(3); break; // 1st T
                    default  : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 17: // Di-nucleotides 2nd position (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 2) {
                switch (toupper(seq[i+1])) {
                    case 'A' : seq_data.push_back(0); break; // 2nd A
                    case 'C' : seq_data.push_back(1); break; // 2nd C
                    case 'G' : seq_data.push_back(2); break; // 2nd G
                    case 'T' : seq_data.push_back(3); break; // 2nd T
                    default  : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 18: // Codons (AAA|AAC|...|TTG|TTT)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case 000: seq_data.push_back(0);  break; // AAA
                    case 001: seq_data.push_back(1);  break; // AAC
                    case 002: seq_data.push_back(2);  break; // AAG
                    case 003: seq_data.push_back(3);  break; // AAT
                    case 010: seq_data.push_back(4);  break; // ACA
                    case 011: seq_data.push_back(5);  break; // ACC
                    case 012: seq_data.push_back(6);  break; // ACG
                    case 013: seq_data.push_back(7);  break; // ACT
                    case 020: seq_data.push_back(8);  break; // AGA
                    case 021: seq_data.push_back(9);  break; // AGC
                    case 022: seq_data.push_back(10); break; // AGG
                    case 023: seq_data.push_back(11); break; // AGT
                    case 030: seq_data.push_back(12); break; // ATA
                    case 031: seq_data.push_back(13); break; // ATC
                    case 032: seq_data.push_back(14); break; // ATG
                    case 033: seq_data.push_back(15); break; // ATT
                    case 100: seq_data.push_back(16); break; // CAA
                    case 101: seq_data.push_back(17); break; // CAC
                    case 102: seq_data.push_back(18); break; // CAG
                    case 103: seq_data.push_back(19); break; // CAT
                    case 110: seq_data.push_back(20); break; // CCA
                    case 111: seq_data.push_back(21); break; // CCC
                    case 112: seq_data.push_back(22); break; // CCG
                    case 113: seq_data.push_back(23); break; // CCT
                    case 120: seq_data.push_back(24); break; // CGA
                    case 121: seq_data.push_back(25); break; // CGC
                    case 122: seq_data.push_back(26); break; // CGG
                    case 123: seq_data.push_back(27); break; // CGT
                    case 130: seq_data.push_back(28); break; // CTA
                    case 131: seq_data.push_back(29); break; // CTC
                    case 132: seq_data.push_back(30); break; // CTG
                    case 133: seq_data.push_back(31); break; // CTT
                    case 200: seq_data.push_back(32); break; // GAA
                    case 201: seq_data.push_back(33); break; // GAC
                    case 202: seq_data.push_back(34); break; // GAG
                    case 203: seq_data.push_back(35); break; // GAT
                    case 210: seq_data.push_back(36); break; // GCA
                    case 211: seq_data.push_back(37); break; // GCC
                    case 212: seq_data.push_back(38); break; // GCG
                    case 213: seq_data.push_back(39); break; // GCT
                    case 220: seq_data.push_back(40); break; // GGA
                    case 221: seq_data.push_back(41); break; // GGC
                    case 222: seq_data.push_back(42); break; // GGG
                    case 223: seq_data.push_back(43); break; // GGT
                    case 230: seq_data.push_back(44); break; // GTA
                    case 231: seq_data.push_back(45); break; // GTC
                    case 232: seq_data.push_back(46); break; // GTG
                    case 233: seq_data.push_back(47); break; // GTT
                    case 300: seq_data.push_back(48); break; // TAA
                    case 301: seq_data.push_back(49); break; // TAC
                    case 302: seq_data.push_back(50); break; // TAG
                    case 303: seq_data.push_back(51); break; // TAT
                    case 310: seq_data.push_back(52); break; // TCA
                    case 311: seq_data.push_back(53); break; // TCC
                    case 312: seq_data.push_back(54); break; // TCG
                    case 313: seq_data.push_back(55); break; // TCT
                    case 320: seq_data.push_back(56); break; // TGA
                    case 321: seq_data.push_back(57); break; // TGC
                    case 322: seq_data.push_back(58); break; // TGG
                    case 323: seq_data.push_back(59); break; // TGT
                    case 330: seq_data.push_back(60); break; // TTA
                    case 331: seq_data.push_back(61); break; // TTC
                    case 332: seq_data.push_back(62); break; // TTG
                    case 333: seq_data.push_back(63); break; // TTT
                    default:  seq_data.push_back(64); break; // In case of other characters
                }
            }
            break;
        case 19: // Codons 1st + 2nd positions (AA|AC|...|TG|TT)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(2,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case 00: seq_data.push_back(0);  break; // AA
                    case 01: seq_data.push_back(1);  break; // AC
                    case 02: seq_data.push_back(2);  break; // AG
                    case 03: seq_data.push_back(3);  break; // AT
                    case 10: seq_data.push_back(4);  break; // CA
                    case 11: seq_data.push_back(5);  break; // CC
                    case 12: seq_data.push_back(6);  break; // CG
                    case 13: seq_data.push_back(7);  break; // CT
                    case 20: seq_data.push_back(8);  break; // GA
                    case 21: seq_data.push_back(9);  break; // GC
                    case 22: seq_data.push_back(10); break; // GG
                    case 23: seq_data.push_back(11); break; // GT
                    case 30: seq_data.push_back(12); break; // TA
                    case 31: seq_data.push_back(13); break; // TC
                    case 32: seq_data.push_back(14); break; // TG
                    case 33: seq_data.push_back(15); break; // TT
                    default: seq_data.push_back(16); break; // In case of other characters
                }
            }
            break;
        case 20: // Codons 1st + 3rd positions (AA|AC|...|TG|TT)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(1,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case 00: seq_data.push_back(0);  break; // AA
                    case 01: seq_data.push_back(1);  break; // AC
                    case 02: seq_data.push_back(2);  break; // AG
                    case 03: seq_data.push_back(3);  break; // AT
                    case 10: seq_data.push_back(4);  break; // CA
                    case 11: seq_data.push_back(5);  break; // CC
                    case 12: seq_data.push_back(6);  break; // CG
                    case 13: seq_data.push_back(7);  break; // CT
                    case 20: seq_data.push_back(8);  break; // GA
                    case 21: seq_data.push_back(9);  break; // GC
                    case 22: seq_data.push_back(10); break; // GG
                    case 23: seq_data.push_back(11); break; // GT
                    case 30: seq_data.push_back(12); break; // TA
                    case 31: seq_data.push_back(13); break; // TC
                    case 32: seq_data.push_back(14); break; // TG
                    case 33: seq_data.push_back(15); break; // TT
                    default: seq_data.push_back(16); break; // In case of other characters
                }
            }
            break;
        case 21: // Codons 2nd + 3rd positions (AA|AC|...|TG|TT)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(0,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case 00: seq_data.push_back(0);  break; // AA
                    case 01: seq_data.push_back(1);  break; // AC
                    case 02: seq_data.push_back(2);  break; // AG
                    case 03: seq_data.push_back(3);  break; // AT
                    case 10: seq_data.push_back(4);  break; // CA
                    case 11: seq_data.push_back(5);  break; // CC
                    case 12: seq_data.push_back(6);  break; // CG
                    case 13: seq_data.push_back(7);  break; // CT
                    case 20: seq_data.push_back(8);  break; // GA
                    case 21: seq_data.push_back(9);  break; // GC
                    case 22: seq_data.push_back(10); break; // GG
                    case 23: seq_data.push_back(11); break; // GT
                    case 30: seq_data.push_back(12); break; // TA
                    case 31: seq_data.push_back(13); break; // TC
                    case 32: seq_data.push_back(14); break; // TG
                    case 33: seq_data.push_back(15); break; // TT
                    default: seq_data.push_back(16); break; // In case of other characters
                }
            }
            break;
        case 22: // Codons 1st + 2nd positions (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    if (j == i || j == i + 1) {
                        switch (toupper(seq[j])) {
                            case 'A': seq_data.push_back(0); break;
                            case 'C': seq_data.push_back(1); break;
                            case 'G': seq_data.push_back(2); break;
                            case 'T': seq_data.push_back(3); break;
                            case 'U': seq_data.push_back(3); break;
                            default : seq_data.push_back(4); break;
                        }
                    }
                }
            }
            break;
        case 23: // Codons 1st + 3rd positions (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    if (j == i || j == i + 2) {
                        switch (toupper(seq[j])) {
                            case 'A': seq_data.push_back(0); break;
                            case 'C': seq_data.push_back(1); break;
                            case 'G': seq_data.push_back(2); break;
                            case 'T': seq_data.push_back(3); break;
                            case 'U': seq_data.push_back(3); break;
                            default : seq_data.push_back(4); break;
                        }
                    }
                }
            }
            break;
        case 24: // Codons 2nd + 3rd positions (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    if (j == i + 1 || j == i + 2) {
                        switch (toupper(seq[j])) {
                            case 'A': seq_data.push_back(0); break;
                            case 'C': seq_data.push_back(1); break;
                            case 'G': seq_data.push_back(2); break;
                            case 'T': seq_data.push_back(3); break;
                            case 'U': seq_data.push_back(3); break;
                            default : seq_data.push_back(4); break;
                        }
                    }
                }
            }
            break;
        case 25: // Codons 1st position (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(2,1);
                triplet.erase(1,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case   0: seq_data.push_back(0); break;
                    case   1: seq_data.push_back(1); break;
                    case   2: seq_data.push_back(2); break;
                    case   3: seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 26: // Codons 2nd position (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(2,1);
                triplet.erase(0,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case   0: seq_data.push_back(0); break;
                    case   1: seq_data.push_back(1); break;
                    case   2: seq_data.push_back(2); break;
                    case   3: seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 27: // Codons 3rd position (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(1,1);
                triplet.erase(0,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case   0: seq_data.push_back(0); break;
                    case   1: seq_data.push_back(1); break;
                    case   2: seq_data.push_back(2); break;
                    case   3: seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 28: // 10-state genotype data
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    case 'K': seq_data.push_back(4); break;
                    case 'M': seq_data.push_back(5); break;
                    case 'R': seq_data.push_back(6); break;
                    case 'Y': seq_data.push_back(7); break;
                    case 'S': seq_data.push_back(8); break;
                    case 'W': seq_data.push_back(9); break;
                    default : seq_data.push_back(10);break; // In case of other characters
                }
            }
            break;
        case 29: // 14-state genotype data
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    case 'K': seq_data.push_back(4); break;
                    case 'M': seq_data.push_back(5); break;
                    case 'R': seq_data.push_back(6); break;
                    case 'Y': seq_data.push_back(7); break;
                    case 'S': seq_data.push_back(8); break;
                    case 'W': seq_data.push_back(9); break;
                    case 'B': seq_data.push_back(10);break;
                    case 'D': seq_data.push_back(11);break;
                    case 'H': seq_data.push_back(12);break;
                    case 'V': seq_data.push_back(13);break;
                    default : seq_data.push_back(14);break; // In case of other characters
                }
            }
            break;
        case 30: // amino acids (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'D': seq_data.push_back(2); break;
                    case 'E': seq_data.push_back(3); break;
                    case 'F': seq_data.push_back(4); break;
                    case 'G': seq_data.push_back(5); break;
                    case 'H': seq_data.push_back(6); break;
                    case 'I': seq_data.push_back(7); break;
                    case 'K': seq_data.push_back(8); break;
                    case 'L': seq_data.push_back(9); break;
                    case 'M': seq_data.push_back(10);break;
                    case 'N': seq_data.push_back(11);break;
                    case 'P': seq_data.push_back(12);break;
                    case 'Q': seq_data.push_back(13);break;
                    case 'R': seq_data.push_back(14);break;
                    case 'S': seq_data.push_back(15);break;
                    case 'T': seq_data.push_back(16);break;
                    case 'V': seq_data.push_back(17);break;
                    case 'W': seq_data.push_back(18);break;
                    case 'Y': seq_data.push_back(19);break;
                    default : seq_data.push_back(20);break; // In case of other characters
                }
            }
            break;
        default: // Dayhoff-6 (AGPST|DENQ|HKR|MIVL|WFY|C)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(0); break;
                    case 'P': seq_data.push_back(0); break;
                    case 'S': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(0); break;
                    case 'D': seq_data.push_back(1); break;
                    case 'E': seq_data.push_back(1); break;
                    case 'N': seq_data.push_back(1); break;
                    case 'Q': seq_data.push_back(1); break;
                    case 'H': seq_data.push_back(2); break;
                    case 'K': seq_data.push_back(2); break;
                    case 'R': seq_data.push_back(2); break;
                    case 'M': seq_data.push_back(3); break;
                    case 'I': seq_data.push_back(3); break;
                    case 'L': seq_data.push_back(3); break;
                    case 'V': seq_data.push_back(3); break;
                    case 'F': seq_data.push_back(4); break;
                    case 'W': seq_data.push_back(4); break;
                    case 'Y': seq_data.push_back(4); break;
                    case 'C': seq_data.push_back(5); break;
                    default : seq_data.push_back(6); break; // In case of other characters
                }
            }
            break;
    }
    return(seq_data);
}


// Function that reads input file and stores data in two 2D containers
void Read_Input(string inname, unsigned datatype){
    unsigned long alignment_length(0);
    unsigned long counter(0);
    string seq(""), str(""), tmp(""); // temporary string used to store input
    vector<int> sequence;     // temporary vector used to store input
    ifstream infile;
    
    infile.open(inname.c_str());
    if (!infile) {
        cerr << "Program aborted due to error: input file not found" << endl;
        exit(1);
    }
    while (getline(infile, str)) {
        if (!str.empty()) {
            // remove blank space in string
            tmp.clear();
            for (std::string::size_type i = 0; i != str.size(); ++i) {
                if (!isblank(str[i])) {
                    tmp.push_back(str[i]);
                }
            }
            if (tmp[0] == '>') {
                if (seq.size() > 0) {
                    if (datatype > 14 && datatype < 18) {
                        if (seq.size() % 2 != 0) {
                            std::cerr << "\nERROR: expected sequence of di-nucleotides" << "\n" << std::endl;
                            exit(1);
                        }
                    }
                    if (datatype > 17 && datatype < 28) {
                        if (seq.size() % 3 != 0) {
                            std::cerr << "\nERROR: expected sequence of codons" << "\n" << std::endl;
                            exit(1);
                        }
                    }
                    sequence = Translator(datatype, seq);
                    alignment.push_back(sequence); // stores sequence in vector
                    if (alignment_length == 0)
                        alignment_length = sequence.size();
                    sequence.clear();
                    seq.clear();
                }
                tmp.erase(tmp.begin()); // removes first character from name
                taxon.push_back(tmp); // stores sequence name in vector
            } else {
                seq += tmp;
            }
            str.clear();
        }
    }
    // Store last sequence in vector
    if (seq.size() > 0) {
        if (datatype > 14 && datatype < 18) {
            if (seq.size() % 2 != 0) {
                cerr << "\nProgram aborted due to error: expected sequence of di-nucleotides" << "\n" << endl;
                exit(1);
            }
        }
        if (datatype > 17 && datatype < 28) {
            if (seq.size() % 3 != 0) {
                cerr << "\nProgram aborted due to error: expected sequence of codons" << "\n" << endl;
                exit(1);
            }
        }
        sequence = Translator(datatype, seq);
        alignment.push_back(sequence);
    } else {
        cerr << "Program aborted due to error: last sequence empty" << "\n" << endl;
        exit(1);
    }
    //Check whether the sequence names are unique
    for (vector<string>::const_iterator iter1 = taxon.begin(); iter1 != taxon.end(); ++iter1) {
        for (vector<string>::const_iterator iter2 = iter1 + 1; iter2 != taxon.end(); ++iter2) {
            if (*iter1 == *iter2) {
                cerr << "Program aborted due to error: Sequence name not unique -- look for " << *iter1 << "\n" << endl;
                exit(1);
            }
        }
    }
    // Check whether the sequences have the same length
    for (vector<vector<int> >::const_iterator iter = alignment.begin()+1; iter != alignment.end(); ++iter) {
        ++counter;
        sequence = *iter;
        if (sequence.size() != alignment_length) {
            cerr << "Program aborted due to error: sequences 1 and " << counter << " differ in length!\n" << endl;
            exit(1);
        }
    }
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// NOTE PERTAINING TO THE FOLLOWING FOUR FUNCTIONS                            //
//                                                                            //
// Source: gaussian_distribution_tail.c                                       //
// Source: chi-square_distribution_tail.c                                     //
// Author: Dick Horn (mathretprogr@gmail.com)                                 //
// Note: Used with permission from the author (Wednesday, 9 July 2014 4:30)   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


// This function returns the probability that a random variable with a standard
// Normal (Gaussian) distribution has a value greater than "x"
long double xGaussian_Distribution_Tail( long double x ) {
    long double sqrt2 = 0.7071067811865475244008443621048490L;
    return  0.5L * erfcl(sqrt2 * x );
}

// The number of degrees of freedom, nu, is an even integer, nu = 2*n.
// The argument x is chi^2 / 2.
static long double Sum_Poisson_Terms(long double x, int n) {
    int k;
    long double term;
    long double sum;
    
    term = 1.0L;
    sum = 1.0L;
    for (k = 1; k < n; k++) {
        term *= (x / k);
        sum += term;
    }
    return expl(-x)*sum;
}

static long double Sum_Over_Odd_Terms(long double x, int dof) {
    int k;
    int n;
    long double term;
    long double sum;
    long double sqrtx;
    long double twooverpi;
    
    twooverpi = 0.6366197723675813430755350534900574L;
    
    if (dof == 1) return 2.0L * xGaussian_Distribution_Tail( sqrtl(x) );
    n = (dof - 1) / 2;
    sqrtx = sqrtl(x);
    term = sqrtx;
    sum = sqrtx;
    for (k = 2; k <=n; k++) {
        term *= ( x / (k + k - 1) );
        sum += term;
    };
    return 2.0L * xGaussian_Distribution_Tail(sqrtx) + sqrtl(twooverpi) * expl(-x/2.0L)*sum;
}

// This function returns the probability that a random chi-squared-distributed
// variable has a value greater than "x"
long double xChi_Square_Distribution_Tail(long double x, int dof) {
    
    if (dof <= 0) return 0.0L;
    
    if (dof % 2 == 0)
        return Sum_Poisson_Terms(x/2.0L, dof/2);
    else
        return Sum_Over_Odd_Terms(x,dof);
}


//MAIN PROGRAM THAT CALLS THE VARIOUS FUNCTIONS ABOVE
int main(int argc, char** argv){
    unsigned rows_columns(0);
    unsigned dataType(0);          // Datatype determines the size of matrices used
    unsigned df(0);                // Degrees of freedom
    unsigned long dm[max_array][max_array]; // 2D divergence matrix
    unsigned long sum_dm(0);       // Sum of elements in divergence matrix
    unsigned long sum_diag_dm(0);  // Sum of diagonal elements in divergence matrix
    unsigned long rejectBonf(0);   // Rejected tests (Bonferroni 1936)
    unsigned long rejectHolm(0);   // Rejected tests (Holm 1979)
    unsigned long rejectBenYek(0); // Rejected tests (Benjamini & Yekutieli 2001)
    unsigned long below_Tau(0);    // P values below tau
    unsigned long counter(0), total;
    double threshold(0.0);     // Threshold used in statistical tests
    double adjustment(0.0);    // Variable used in Benjamini & Yekutieli's (2001) test
    double p_Distance;         // p-distance
    double d_EuclideanFS(0.0); // Euclidean distance -- full symmetry
    double d_EuclideanMS(0.0); // Euclidean distance -- marginal symmetry
    double d_Bowker(0.0);      // Compositional distance -- full symmetry
    double row_sum[max_array], col_sum[max_array];
    double min_Prob(1.0);      //, family_Wise_Error_Rate(0.05);
    double max_p_Dist(-numeric_limits<double>::max());
    double min_p_Dist(numeric_limits<double>::max());
    double max_delta_B(-numeric_limits<double>::max());
    double min_delta_B(numeric_limits<double>::max());
    double max_delta_EuclFS(-numeric_limits<double>::max());
    double min_delta_EuclFS(numeric_limits<double>::max());
    double max_delta_EuclMS(-numeric_limits<double>::max());
    double min_delta_EuclMS(numeric_limits<double>::max());
    long double BMPTS(0.0);   // Bowker's (1948) test statistic, initialised
    long double P_value(1.0); // The observed P value, initialised
    double Bonferroni;        // Family-wise error rate (Bonferroni)
    vector<int> length;       // Length of alignment
    vector<int> row_of_int;
    vector<unsigned> row_of_unsigned;
    vector<unsigned long> row_of_unsigned_long;
    vector<double> row_of_double;
    vector<size_t> seq1;              // Vector holding sequence 1 identifier
    vector<size_t> seq2;              // Vector holding sequence 2 identifier
    vector<unsigned long> sites;      // Vector holding the number of sites compared
    vector<double> test_Bowker;       // Vector holding the test statistic from Bowker (1948)
    vector<unsigned> degrees_Freedom; // Vector holding the degrees of freedom
    vector<double> P_Bowker;          // Vector holding P values from Bowker's (1948) test
    vector<double> P_exp;             // Vector holding expected P values
    vector<size_t> rank;              // Vector holding the rank according to P value (low - high)
    vector<size_t> order_P;           // Vector holding the order in which sequences were compared
    vector<double> p_Dist;            // Vector holding the p-distance
    vector<double> delta_B;           // Vector holding compositional distance based on Bowker (1948)
    vector<double> d_EuclidFS;        // Vector holding Euclidean distance based on full symmetry)
    vector<double> d_EuclidMS;        // Vector holding Euclidean distance based on marginal symmetry)
    vector<double> Bonf;              // Vector holding family-wise error rate (Bonferroni 1936)
    vector<double> Holm;              // Vector holding family-wise error rate (Holm 1979)
    vector<double> BenYek;            // Vector holding false discovery rate (Benjamini & Yekulieti 2001)
    vector<vector<double> > mat_Bowker;  // Matrix initialised with 0.0
    vector<vector<double> > mat_Prob;    // Matrix holding probabilities
    vector<vector<double> > mat_pDist;   // Matrix holding p-distances
    vector<vector<double> > mat_delta_B; // Matrix holding compositional distance (Bowker)
    vector<vector<double> > mat_defs;    // Matrix holding compositional distance (Euklidean; full symmetry)
    vector<vector<double> > mat_dems;    // Matrix holding compositional distance (Euklidean; marginal symmetry)
    vector<vector<unsigned> > mat_Bonferroni; // Matrix holding 0 or 1 (fail or accept)
    vector<vector<unsigned> > mat_Holm;       // Matrix holding 0 or 1 (fail or accept)
    vector<vector<unsigned> > mat_BenYek;     // Matrix holding 0 or 1 (fail or accept)
    vector<vector<unsigned> > mat_df;    // Matrix holding degrees of freedom
    vector<vector<unsigned long> > mat_sites; // Matrix holding number of sites
    string name_1, name_2; // temporary variables holding names of two sequences
    string survey, str;
    string inName, outName1, outName2, outName3, outName4, outName5, outName6, outName7, outName8, outName9;
    ofstream outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9;
    
    if(argc != 5) {
        cerr << "\nHomo v" << VERSION << " Copyright 2019-22, Lars Jermiin" << endl;
        cerr << " Contact: lars.jermiin [at] anu.edu.au / ucd.ie" << endl;
        cerr << "\nERROR -- use command: homo <infile> <b|f> <1|...|31> <error rate>\n" << endl;
        cerr << "  infile   Fasta-formatted alignment" << endl;
        cerr << "     b|f   Brief or full report of results" << endl;
        cerr << "       1   Nucleotides; 4 states (A|C|G|T)" << endl;
        cerr << "       2   Nucleotides; 3 states (C|T|AG)" << endl;
        cerr << "       3   Nucleotides; 3 states (A|G|CT)" << endl;
        cerr << "       4   Nucleotides; 3 states (A|T|CG)" << endl;
        cerr << "       5   Nucleotides; 3 states (C|G|AT)" << endl;
        cerr << "       6   Nucleotides; 3 states (A|C|GT)" << endl;
        cerr << "       7   Nucleotides; 3 states (G|T|AC)" << endl;
        cerr << "       8   Nucleotides; 2 states (GT|AC)" << endl;
        cerr << "       9   Nucleotides; 2 states (AG|CT)" << endl;
        cerr << "      10   Nucleotides; 2 states (GC|AT)" << endl;
        cerr << "      11   Nucleotides; 2 states (A|CGT)" << endl;
        cerr << "      12   Nucleotides; 2 states (C|AGT)" << endl;
        cerr << "      13   Nucleotides; 2 states (G|ACT)" << endl;
        cerr << "      14   Nucleotides; 2 states (T|ACG)" << endl;
        cerr << "      15   Di-nucleotides (pos 12); 16 states (AA|AC|...|TG|TT)" << endl;
        cerr << "      16   Di-nucleotides (pos  1);  4 states (A|C|G|T)" << endl;
        cerr << "      17   Di-nucleotides (pos  2);  4 states (A|C|G|T)" << endl;
        cerr << "      18   Codons (pos 123); 64 states (AAA|AAC|...|TTG|TTT)" << endl;
        cerr << "      19   Codons (pos  12); 16 states (AA|AC|...|TG|TT)" << endl;
        cerr << "      20   Codons (pos  13); 16 states (AA|AC|...|TG|TT)" << endl;
        cerr << "      21   Codons (pos  23); 16 states (AA|AC|...|TG|TT)" << endl;
        cerr << "      22   Codons (pos  12);  4 states (A|C|G|T)" << endl;
        cerr << "      23   Codons (pos  13);  4 states (A|C|G|T)" << endl;
        cerr << "      24   Codons (pos  23);  4 states (A|C|G|T)" << endl;
        cerr << "      25   Codons (pos   1);  4 states (A|C|G|T)" << endl;
        cerr << "      26   Codons (pos   2);  4 states (A|C|G|T)" << endl;
        cerr << "      27   Codons (pos   3);  4 states (A|C|G|T)" << endl;
        cerr << "      28   Genotypes; 10 states (A|C|G|T|K|M|R|Y|S|W)" << endl;
        cerr << "      29   Genotypes; 14 states (A|C|G|T|K|M|R|Y|S|W|B|D|H|V)" << endl;
        cerr << "      30   Amino acids; 20 states (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C)" << endl;
        cerr << "      31   Amino acids;  6 states (AGPST|DENQ|HKR|MIVL|WFY|C) [D6]" << endl;
        cerr << endl;
        exit(1);
    }
    inName = argv[1];
    survey = argv[2];
    dataType = stoi(argv[3]);
    threshold = stod(argv[4]);
    if (dataType < 1 || dataType > 31) {
        cerr << "\nPROGRAM ABORTED - incorrect choice of data: [1|...|31]\n" << endl;
        exit(1);
    }
    if (toupper(survey[0]) != 'F' && toupper(survey[0]) != 'B') {
        cerr << "\nPROGRAM ABORTED - incorrect choice of output: [b|f]\n" << endl;
        exit(1);
    }
    if (threshold < 0.0 || threshold > 1.0) {
        cerr << "\nPROGRAM ABORTED - incorrect choice of error rate: [0.0 - 1.0]\n" << endl;
        exit(1);
    }
    if (toupper(survey[0]) == 'F') {
        outName1.clear();
        for (string::size_type i = 0; i != inName.size() && inName[i] != '.'; ++i) {
            outName1 += inName[i];
        }
        outName9 = outName1 + "_delta_EMS.nex";
        outName8 = outName1 + "_delta_EFS.nex";
        outName7 = outName1 + "_delta_Bowker.nex";
        outName6 = outName1 + "_p-distances.nex";
        outName5 = outName1 + "_Benja-Yekut.csv";
        outName4 = outName1 + "_Holm.csv";
        outName3 = outName1 + "_Bonferroni.csv";
        outName2 = outName1 + "_P-values.csv";
        outName1 = outName1 + "_Summary.csv";
    }
    switch (dataType) {
        case  1: rows_columns = FOUR; break;
        case  2: rows_columns = THREE; break;
        case  3: rows_columns = THREE; break;
        case  4: rows_columns = THREE; break;
        case  5: rows_columns = THREE; break;
        case  6: rows_columns = THREE; break;
        case  7: rows_columns = THREE; break;
        case  8: rows_columns = TWO; break;
        case  9: rows_columns = TWO; break;
        case 10: rows_columns = TWO; break;
        case 11: rows_columns = TWO; break;
        case 12: rows_columns = TWO; break;
        case 13: rows_columns = TWO; break;
        case 14: rows_columns = TWO; break;
        case 15: rows_columns = SIXTEEN; break;
        case 16: rows_columns = FOUR; break;
        case 17: rows_columns = FOUR; break;
        case 18: rows_columns = SIXTYFOUR; break;
        case 19: rows_columns = SIXTEEN; break;
        case 20: rows_columns = SIXTEEN; break;
        case 21: rows_columns = SIXTEEN; break;
        case 22: rows_columns = FOUR; break;
        case 23: rows_columns = FOUR; break;
        case 24: rows_columns = FOUR; break;
        case 25: rows_columns = FOUR; break;
        case 26: rows_columns = FOUR; break;
        case 27: rows_columns = FOUR; break;
        case 28: rows_columns = TEN; break;
        case 29: rows_columns = FOURTEEN; break;
        case 30: rows_columns = TWENTY; break;
        default: rows_columns = SIX; break;
    }
    Read_Input(inName, dataType);
    length = alignment[0];
    // Initialising different types of rows in matrices
    for (vector<int>::size_type i = 0; i != taxon.size(); i++) {
        row_of_unsigned.push_back(0);
        row_of_unsigned_long.push_back(1);
        row_of_double.push_back(0.0);
    }
    for (vector<vector<int> >::size_type i = 0; i != taxon.size(); i++) {
        mat_df.push_back(row_of_unsigned); // Matrix of degrees of freedom
        mat_sites.push_back(row_of_unsigned_long); // Matrix of sites compared
        mat_Bowker.push_back(row_of_double); // Matrix of Bowker's test statistic
    }
    mat_Prob = mat_Bowker; // Matrix of probabilities
    mat_pDist = mat_Bowker; // Matrix of p distances
    mat_delta_B = mat_Bowker; // Matrix of compositional distances (Bowker)
    mat_defs = mat_Bowker;
    mat_dems = mat_Bowker;
    mat_Bonferroni = mat_df;
    mat_Holm = mat_df;
    mat_BenYek = mat_df;
    if (toupper(survey[0]) == 'F') {
        outfile1.open(outName1.c_str());
        outfile1 << "Program, Homo" << endl;
        outfile1 << "Version, " << VERSION << endl;
        outfile1 << "Input file," << inName << endl;
        outfile1 << "Characters,";
        switch (dataType) {
            case  1: outfile1 << "Nucleotides (A|C|G|T)" << endl; break;
            case  2: outfile1 << "Nucleotides recoded (C|T|AG)" << endl; break;
            case  3: outfile1 << "Nucleotides recoded (A|G|CT)" << endl; break;
            case  4: outfile1 << "Nucleotides recoded (A|T|CG)" << endl; break;
            case  5: outfile1 << "Nucleotides recoded (C|G|AT)" << endl; break;
            case  6: outfile1 << "Nucleotides recoded (A|C|GT)" << endl; break;
            case  7: outfile1 << "Nucleotides recoded (G|T|AC)" << endl; break;
            case  8: outfile1 << "Nucleotides recoded (GT|AC)" << endl; break;
            case  9: outfile1 << "Nucleotides recoded (AG|CT)" << endl; break;
            case 10: outfile1 << "Nucleotides recoded (GC|AT)" << endl; break;
            case 11: outfile1 << "Nucleotides recoded (A|CGT)" << endl; break;
            case 12: outfile1 << "Nucleotides recoded (C|AGT)" << endl; break;
            case 13: outfile1 << "Nucleotides recoded (G|ACT)" << endl; break;
            case 14: outfile1 << "Nucleotides recoded (T|ACG)" << endl; break;
            case 15: outfile1 << "Di-nucleotides (AA|AC|...|TG|TT)" << endl; break;
            case 16: outfile1 << "Di-nucleotides 1st position (A|C|G|T)" << endl; break;
            case 17: outfile1 << "Di-nucleotides 2nd position (A|C|G|T)" << endl; break;
            case 18: outfile1 << "Codons (AAA|AAC|...|TTG|TTT)" << endl; break;
            case 19: outfile1 << "Codons 1st + 2nd positions (AA|AC|...|TG|TT)" << endl; break;
            case 20: outfile1 << "Codons 1st + 3rd positions (AA|AC|...|TG|TT)" << endl; break;
            case 21: outfile1 << "Codons 2nd + 3rd positions (AA|AC|...|TG|TT)" << endl; break;
            case 22: outfile1 << "Codons 1st + 2nd positions (A|C|G|T)" << endl; break;
            case 23: outfile1 << "Codons 1st + 3rd positions (A|C|G|T)" << endl; break;
            case 24: outfile1 << "Codons 2nd + 3rd positions (A|C|G|T)" << endl; break;
            case 25: outfile1 << "Codons 1st position (A|C|G|T)" << endl; break;
            case 26: outfile1 << "Codons 2nd position (A|C|G|T)" << endl; break;
            case 27: outfile1 << "Codons 3rd position (A|C|G|T)" << endl; break;
            case 28: outfile1 << "Genotypes (A|C|G|T|K|M|R|Y|S|W)" << endl; break;
            case 29: outfile1 << "Genotypes (A|C|G|T|K|M|R|Y|S|W|B|D|H|V)" << endl; break;
            case 30: outfile1 << "Amino acids (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C)" << endl; break;
            default: outfile1 << "Recoded amino acids (AGPST|DENQ|HKR|MIVL|WFY|C) [D6]" << endl; break;
        }
        outfile1 << "Sequences," << taxon.size() << endl << endl;
        outfile1 << "Taxon1,Taxon2,Sites,Bowker,df,P_obs,Rank,P_exp,P_Bonferroni,Test,P_Holm,Test,P_Benjamini_Yekutieli,Test,d_Bowker,d_Euclidean_FS,d_Euclidean_MS,p_distance" << endl;
        cout << endl;
    }
    total = taxon.size() * (taxon.size() - 1)/2;
    counter = total;

    // Loops generating primary estimates (i.e., P values, p-distances, and compositional distances)
    for (vector<vector<int> >::size_type iter1 = 0; iter1 != alignment.size(); iter1++) {
        for (vector<vector<int> >::size_type iter2 = iter1 + 1; iter2 != alignment.size(); iter2++) {
            // Initialise divergence matrix using 0
            for (size_t i = 0; i < max_array; i++) {
                for (size_t j = 0; j < max_array; j++) {
                    dm[i][j] = 0;
                }
            }
            // Generate divergence matrix
            for (vector<int>::size_type i = 0; i < length.size(); i++) {
                dm[alignment[iter1][i]][alignment[iter2][i]]++;
            }
            // Compute the sum of elements in the divergence matrix
            sum_dm = 0;
            for (size_t m = 0; m < rows_columns; m++) {
                for (size_t n = 0; n < rows_columns; n++) {
                    sum_dm += dm[m][n];
                }
            }
            // Storing the number of sites compared for later use
            mat_sites[iter1][iter2] = sum_dm;
            mat_sites[iter2][iter1] = sum_dm;

            // Compute the sum of diagonal elements in the divergence matrix
            sum_diag_dm = 0;
            for (size_t i = 0; i < rows_columns; i++) {
                sum_diag_dm += dm[i][i];
            }
            // Compute the p-distance
            p_Distance = (double)(sum_dm - sum_diag_dm)/(double)sum_dm;
            // Record the largest and smallest p-distance
            if (p_Distance > max_p_Dist) {
                max_p_Dist = p_Distance;
            };
            if (p_Distance < min_p_Dist) {
                min_p_Dist = p_Distance;
            };
            // Store p-distance in matrix for later use
            mat_pDist[iter1][iter2] = p_Distance;
            mat_pDist[iter2][iter1] = p_Distance;

            // Compute the test statistic and degrees of freedom for Bowker's (1948) test
            BMPTS = 0.0;
            df = 0;
            for (size_t m = 0; m < rows_columns; m++) {
                for (size_t n = m+1; n < rows_columns; n++) {
                    if (dm[m][n] + dm[n][m] > 0) {
                        ++df;
                        BMPTS = BMPTS + ((long double)(SQR(dm[m][n] - dm[n][m])))/(dm[m][n] + dm[n][m]);
                    }
                }
            }
            // Store Bowker's test statistic for later use
            mat_Bowker[iter1][iter2] = BMPTS;
            mat_Bowker[iter2][iter1] = BMPTS;
            // Store the degrees of freedom for later use
            mat_df[iter1][iter2] = df;
            mat_df[iter2][iter1] = df;
            // Compute the probability of Bowker's (1948) test statistic
            if (df == 0) {
                P_value = 1.0;
            } else {
                P_value = xChi_Square_Distribution_Tail(BMPTS, df);
            }
            // Count the number of P values below threshold
            if (P_value < threshold) {
                ++below_Tau;
            }
            // Record the smallest P value and the sequence names that produced it
            if (P_value < min_Prob) {
                min_Prob = P_value;
                name_1 = taxon[iter1];
                name_2 = taxon[iter2];
            }
            // Store probability in matrix for later use
            mat_Prob[iter1][iter2] = P_value;
            mat_Prob[iter2][iter1] = P_value;
            // Create vector of P values for sorting and creating a rank
            P_Bowker.push_back(P_value);
            
            // Compute compositional distance -- 1. calculating delta_B
            if (df == 0) {
                d_Bowker = 0.0;
            } else {
//               d_Bowker = sqrt(BMPTS/df);
//               d_Bowker = d_Bowker/sum_dm; // Per site scaling introduced on 21 Oct 2021
                d_Bowker = sqrt(BMPTS/(df * sum_dm)); // Correction on 5 Mar 2022
            }
            // Store compositional distance in matrix for later use
            mat_delta_B[iter1][iter2] = d_Bowker;
            mat_delta_B[iter2][iter1] = d_Bowker;
           // Record the smallest and largest values of delta_B
            if (d_Bowker > max_delta_B) {
                max_delta_B = d_Bowker;
            }
            if (d_Bowker < min_delta_B) {
                min_delta_B = d_Bowker;
            }
            // Compute compositional distance -- 2. calculating Euclidean distance (full sym.)
            d_EuclideanFS = 0.0;
            for (size_t m = 0; m < rows_columns; m++) {
                for (size_t n = m+1; n < rows_columns; n++) {
                    // dm is defined as unsigned long so this if-else statement is needed
                    if (dm[m][n] > dm[n][m]) {
                        d_EuclideanFS = d_EuclideanFS + SQR(((double)(dm[m][n] - dm[n][m]))/sum_dm);
                    } else {
                        d_EuclideanFS = d_EuclideanFS + SQR(((double)(dm[n][m] - dm[m][n]))/sum_dm);
                    }
                }
            }
            d_EuclideanFS = sqrt(d_EuclideanFS);
            // Enter compositional distance into matrix for later use
            mat_defs[iter1][iter2] = d_EuclideanFS;
            mat_defs[iter2][iter1] = d_EuclideanFS;
            // Record the smallest and largest values of delta_Euclidean_FS
            if (d_EuclideanFS > max_delta_EuclFS) {
                max_delta_EuclFS = d_EuclideanFS;
            }
            if (d_EuclideanFS < min_delta_EuclFS) {
                min_delta_EuclFS = d_EuclideanFS;
            }
            // Compute compositional distance -- 3. calculating Euclidean distance (mar. sym.)
            d_EuclideanMS = 0.0;
            // Initialise marginal elements in divergence matrix to zero
            for (size_t i = 0; i < rows_columns; i++) {
                row_sum[i] = 0.0;
                col_sum[i] = 0.0;
            }
            // Compute vectors of marginal frequencies
            for (size_t i = 0; i < rows_columns; i++) {
                for (size_t j = 0; j < rows_columns; j++) {
                    if (i != j) {
                        row_sum[i] += dm[i][j];
                    }
                }
                row_sum[i] = row_sum[i]/sum_dm;
            }
            for (size_t i = 0; i != rows_columns; i++) {
                for (size_t j = 0; j != rows_columns; j++) {
                    if (i != j) {
                        col_sum[i] += dm[j][i];
                    }
                }
                col_sum[i] = col_sum[i]/sum_dm;
            }
            for (size_t i = 0; i < rows_columns; i++) {
                d_EuclideanMS = d_EuclideanMS + (double)(SQR(row_sum[i] - col_sum[i]));
            }
            d_EuclideanMS = sqrt(d_EuclideanMS);
            // Enter compositional distance into matrix for later use
            mat_dems[iter1][iter2] = d_EuclideanMS;
            mat_dems[iter2][iter1] = d_EuclideanMS;
            // Record the smallest and largest values of delta_Euclidean_MS
            if (d_EuclideanMS > max_delta_EuclMS) {
                max_delta_EuclMS = d_EuclideanMS;
            }
            if (d_EuclideanMS < min_delta_EuclMS) {
                min_delta_EuclMS = d_EuclideanMS;
            }
            // Report to terminal how much computation is left
            cout << "\rNumber of comparisons left = " << --counter;
            fflush(NULL);
        }
    }

    // Prepare vectors to rank P values
    for (size_t i = 0; i < total; i++) {
        order_P.push_back(i+1);
        rank.push_back(i+1);
    }
    // Sorting of order_P according to P_values
    for (size_t i = 0; i < total - 1; i++) {
        for (size_t j = i + 1; j < total; j++) {
            // Swap P_values and order_P values
            if (P_Bowker[j] < P_Bowker[i]) {
                long double temp_P;
                temp_P = P_Bowker[i];
                P_Bowker[i] = P_Bowker[j];
                P_Bowker[j] = temp_P;
                // Swap elements in order_P
                unsigned long pair;
                pair = order_P[i];
                order_P[i] = order_P[j];
                order_P[j] = pair;
            }
        }
    }
    // Sorting of rank values according to order_P values
    for (size_t i = 0; i < total - 1; i++) {
        for (size_t j = i + 1; j < total; j++) {
            // Swap P_values and order_P values
            if (order_P[j] < order_P[i]) {
                unsigned long pair;
                pair = order_P[i];
                order_P[i] = order_P[j];
                order_P[j] = pair;
                // Swap elements in rank
                pair = rank[i];
                rank[i] = rank[j];
                rank[j] = pair;
            }
        }
    }

    // Computing adjusted threshold values for FWER and FDR procedures
    adjustment = 0.0;
    for (size_t j = 0; j < total; j++) {
        adjustment = adjustment + 1.0/(j + 1);
    }
    Bonferroni = threshold/total;
    for (size_t i = 0; i < total; i++) {
        P_exp.push_back((double)(i+1)/(total + 1));
        Bonf.push_back(Bonferroni);
        Holm.push_back((double)threshold/(total + 1 - (i+1)));
        BenYek.push_back((double)(threshold * (i+1))/(total * adjustment));
    }

    // Counting rejections
    size_t i(0);
    double cutoffBonf(0.0), cutoffHolm(0.0), cutoffBenYek(0.0);
    for (vector<vector<int> >::size_type iter1 = 0; iter1 != alignment.size(); iter1++) {
        for (vector<vector<int> >::size_type iter2 = iter1 + 1; iter2 != alignment.size(); iter2++) {
            if (iter1 < iter2) {
                if (mat_Prob[iter1][iter2] < Bonferroni) {
                    ++rejectBonf;
                    if (cutoffBonf < Bonferroni) {
                        cutoffBonf = Bonferroni;
                    }
                }
                if (mat_Prob[iter1][iter2] < Holm[rank[i]-1]) {
                    ++rejectHolm;
                    if (cutoffHolm < Holm[rank[i]-1]) {
                        cutoffHolm = Holm[rank[i]-1];
                    }
                }
                if (mat_Prob[iter1][iter2] < BenYek[rank[i]-1]) {
                    ++rejectBenYek;
                    if (cutoffBenYek < BenYek[rank[i]-1]) {
                        cutoffBenYek = BenYek[rank[i]-1];
                    }
                }
                i++;
            }
        }
    }
    // Printing output to terminal and files
    if (toupper(survey[0]) == 'B') {
        // Printing brief summary to terminal
        cout << "\n\nFile,Sequences,Sites,threshold,Tests,P_min,Bonferroni,Holm,Benjamini_Yekutieli" << endl;
        cout << inName << ",";
        cout << taxon.size() << ",";
        cout << length.size() << ",";
        cout << threshold << ",";
        cout << total << ",";
        cout << scientific << min_Prob << ",";
        cout << fixed << rejectBonf << ",";
        cout << fixed << rejectHolm << ",";
        cout << fixed << rejectBenYek << endl;
//        cout << "Offenders " << name_1 << " vs " << name_2 << std::endl;
    } else {
        size_t i(0);
        // Printing full summary to files
        for (vector<vector<int> >::size_type iter1 = 0; iter1 != alignment.size(); iter1++) {
            for (vector<vector<int> >::size_type iter2 = iter1 + 1; iter2 != alignment.size(); iter2++) {
                if (iter1 < iter2) {
                    outfile1 << taxon[iter1] << ",";
                    outfile1 << taxon[iter2] << ",";
                    outfile1 << mat_sites[iter1][iter2] << ",";
                    outfile1 << mat_Bowker[iter1][iter2] << ",";
                    outfile1 << mat_df[iter1][iter2] << ",";
                    outfile1 << mat_Prob[iter1][iter2] << ",";
                    outfile1 << rank[i] << ",";
                    outfile1 << P_exp[rank[i]-1] << ",";
                    outfile1 << Bonferroni << ",";
                    if (mat_Prob[iter1][iter2] < Bonferroni) {
                        outfile1 << "Reject,";
                    } else {
                        outfile1 << "Accept,";
                    }
                    outfile1 << Holm[rank[i]-1] << ",";
                    if (mat_Prob[iter1][iter2] < Holm[rank[i]-1]) {
                        outfile1 << "Reject,";
                    } else {
                        outfile1 << "Accept,";
                    }
                    outfile1 << BenYek[rank[i]-1] << ",";
                    if (mat_Prob[iter1][iter2] < BenYek[rank[i]-1]) {
                        outfile1 << "Reject,";
                    } else {
                        outfile1 << "Accept,";
                    }
                    outfile1 << mat_delta_B[iter1][iter2] << ",";
                    outfile1 << mat_defs[iter1][iter2] << ",";
                    outfile1 << mat_dems[iter1][iter2] << ",";
                    outfile1 << mat_pDist[iter1][iter2] << endl;
                    i++;
                }
            }
        }
        outfile1.close();
        // Printing output intended for .csv files
        outfile2.open(outName2.c_str());
        outfile2 << taxon.size() << endl;
        outfile3.open(outName3.c_str());
        outfile3 << taxon.size() << endl;
        outfile4.open(outName4.c_str());
        outfile4 << taxon.size() << endl;
        outfile5.open(outName5.c_str());
        outfile5 << taxon.size() << endl;
        for (vector<vector<int> >::size_type iter1 = 0; iter1 != alignment.size(); iter1++) {
            outfile2 << taxon[iter1]; // P values
            outfile3 << left << setw(10) << taxon[iter1]; // Bonferroni
            outfile4 << left << setw(10) << taxon[iter1]; // Holm
            outfile5 << left << setw(10) << taxon[iter1]; // Benjamini-Yekutieli
            for (vector<vector<int> >::size_type iter2 = 0; iter2 != alignment.size(); iter2++) {
                outfile2 << "," << fixed << setprecision(11) << mat_Prob[iter1][iter2];
                if (mat_Prob[iter1][iter2] < cutoffBonf) {
                    outfile3 << ",0.0";
                } else {
                    outfile3 << ",1.0";
                }
                if (mat_Prob[iter1][iter2] < cutoffHolm) {
                    outfile4 << ",0.0";
                } else {
                    outfile4 << ",1.0";
                }
                if (mat_Prob[iter1][iter2] < cutoffBenYek) {
                    outfile5 << ",0.0";
                } else {
                    outfile5 << ",1.0";
                }
                i++;
            }
            outfile2 << endl;
            outfile3 << endl;
            outfile4 << endl;
            outfile5 << endl;
        }
        outfile2.close();
        outfile3.close();
        outfile4.close();
        outfile5.close();
        // Printing output intended for .nex files
        outfile6.open(outName6.c_str());
        outfile7.open(outName7.c_str());
        outfile8.open(outName8.c_str());
        outfile9.open(outName9.c_str());
        outfile6 << "#nexus\n" << endl;
        outfile7 << "#nexus\n" << endl;
        outfile8 << "#nexus\n" << endl;
        outfile9 << "#nexus\n" << endl;
        outfile6 << "BEGIN Taxa;" << endl;
        outfile7 << "BEGIN Taxa;" << endl;
        outfile8 << "BEGIN Taxa;" << endl;
        outfile9 << "BEGIN Taxa;" << endl;
        outfile6 << "DIMENSIONS ntax=" << taxon.size() << ";" << endl;
        outfile7 << "DIMENSIONS ntax=" << taxon.size() << ";" << endl;
        outfile8 << "DIMENSIONS ntax=" << taxon.size() << ";" << endl;
        outfile9 << "DIMENSIONS ntax=" << taxon.size() << ";" << endl;
        outfile6 << "TAXLABELS" << endl;
        outfile7 << "TAXLABELS" << endl;
        outfile8 << "TAXLABELS" << endl;
        outfile9 << "TAXLABELS" << endl;
        for (size_t i = 0; i < taxon.size(); i++) {
            outfile6 << "[" << i+1 << "] '" << taxon[i] << "'" << endl;
            outfile7 << "[" << i+1 << "] '" << taxon[i] << "'" << endl;
            outfile8 << "[" << i+1 << "] '" << taxon[i] << "'" << endl;
            outfile9 << "[" << i+1 << "] '" << taxon[i] << "'" << endl;
       }
        outfile6 << ";" << endl;
        outfile7 << ";" << endl;
        outfile8 << ";" << endl;
        outfile9 << ";" << endl;
        outfile6 << "END; [Taxa]\n" << endl;
        outfile7 << "END; [Taxa]\n" << endl;
        outfile8 << "END; [Taxa]\n" << endl;
        outfile9 << "END; [Taxa]\n" << endl;
        outfile6 << "BEGIN Distances;" << endl;
        outfile7 << "BEGIN Distances;" << endl;
        outfile8 << "BEGIN Distances;" << endl;
        outfile9 << "BEGIN Distances;" << endl;
        outfile6 << "DIMENSIONS ntax=" << taxon.size() << ";" << endl;
        outfile7 << "DIMENSIONS ntax=" << taxon.size() << ";" << endl;
        outfile8 << "DIMENSIONS ntax=" << taxon.size() << ";" << endl;
        outfile9 << "DIMENSIONS ntax=" << taxon.size() << ";" << endl;
        outfile6 << "FORMAT labels=no diagonal triangle=both;" << endl;
        outfile7 << "FORMAT labels=no diagonal triangle=both;" << endl;
        outfile8 << "FORMAT labels=no diagonal triangle=both;" << endl;
        outfile9 << "FORMAT labels=no diagonal triangle=both;" << endl;
        outfile6 << "MATRIX" << endl;
        outfile7 << "MATRIX" << endl;
        outfile8 << "MATRIX" << endl;
        outfile9 << "MATRIX" << endl;
        for (vector<vector<int> >::size_type iter1 = 0; iter1 != alignment.size(); iter1++) {
            for (vector<vector<int> >::size_type iter2 = 0; iter2 != alignment.size(); iter2++) {
                outfile6 << " " << fixed << setprecision(11) << mat_pDist[iter1][iter2];
                outfile7 << " " << fixed << setprecision(11) << mat_delta_B[iter1][iter2];
                outfile8 << " " << fixed << setprecision(11) << mat_defs[iter1][iter2];
                outfile9 << " " << fixed << setprecision(11) << mat_dems[iter1][iter2];
            }
            outfile6 << endl;
            outfile7 << endl;
            outfile8 << endl;
            outfile9 << endl;
        }
        outfile6 << ";" << endl;
        outfile7 << ";" << endl;
        outfile8 << ";" << endl;
        outfile9 << ";" << endl;
        outfile6 << "END; [Distances]\n" << endl;
        outfile7 << "END; [Distances]\n" << endl;
        outfile8 << "END; [Distances]\n" << endl;
        outfile9 << "END; [Distances]\n" << endl;
        outfile6.close();
        outfile7.close();
        outfile8.close();
        outfile9.close();
        cout << endl;
        cout << endl;
        cout << "--------------------------------------------------------------------" << endl;
        cout << "   SUMMARY OF ANALYSIS WITH HOMO 2.1" << endl;
        cout << endl;
        cout << "   All estimates ................................ " << outName1 << endl;
        cout << "   Estimates of P values (Bowker 1948) .......... " << outName2 << endl;
        cout << "   Test results (Bonferroni 1936) ............... " << outName3 << endl;
        cout << "   Test results (Holm 1979) ..................... " << outName4 << endl;
        cout << "   Test results (Benjamini & Yekutieli 2001) .... " << outName5 << endl;
        cout << "   p distances .................................. " << outName6 << endl;
        cout << "   Estimates of d_Bowker ........................ " << outName7 << endl;
        cout << "   Estimates of d_Euclidean (Full Sym.) ......... " << outName8 << endl;
        cout << "   Estimates of d_Euclidean (Mar. Sym.) ......... " << outName9 << endl;
        cout << endl;
        cout << "   Positions in alignment ....................... " << length.size() << endl;
        cout << "   Number of tests .............................. " << total << endl;
        cout << "   Smallest P value ............................. " << scientific << min_Prob << endl;
        cout << "   Level of significance (tau) .................. " << fixed << (double) threshold << endl;
        cout << "   Proportion of P values below tau ............. " << fixed << (double)below_Tau/total << endl;
        cout << "   Tests rejected (Bonferroni 1936) ............. " << fixed << rejectBonf << endl;
        cout << "   Tests rejected (Holm 1979) ................... " << fixed << rejectHolm << endl;
        cout << "   Tests rejected (Benjamini & Yekutieli 2001) .. " << fixed << rejectBenYek << endl;
        cout << "   Min(delta_Bowker) ............................ " << fixed << min_delta_B << endl;
        cout << "   Max(delta_Bowker) ............................ " << fixed << max_delta_B << endl;
        cout << "   Min(delta_EuclideanFS) ....................... " << fixed << min_delta_EuclFS << endl;
        cout << "   Max(delta_EuclideanFS) ....................... " << fixed << max_delta_EuclFS << endl;
        cout << "   Min(delta_EuclideanMS) ....................... " << fixed << min_delta_EuclMS << endl;
        cout << "   Max(delta_EuclideanMS) ....................... " << fixed << max_delta_EuclMS << endl;
        cout << "   Min(p distance) .............................. " << fixed << min_p_Dist << endl;
        cout << "   Max(p distance) .............................. " << fixed << max_p_Dist << endl;
        
        if (rejectBonf > 0 || rejectHolm > 0 || rejectBenYek > 0) {
            cout << endl;
            cout << "WARNING:" << endl << endl;
            cout << "   At least one pair of sequences is unlikely to have evolved" << endl;
            cout << "   under the same Markovian process. For further details, see" << endl;
            cout << "   " << outName1 << endl;
        }
        cout << "--------------------------------------------------------------------" << endl;
        cout << endl;
    }
    return 0;
}
