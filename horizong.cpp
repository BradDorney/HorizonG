/*  
 * Author: Brad Dorney
 *
 * Based on NCBI BLAST+ code and API from the NCBI C++ Toolkit
 * http://www.ncbi.nlm.nih.gov/toolkit
 *
 * File Description:
 *   Application for running a BLAST search specialized for HGT analysis.
 *
 */

#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <corelib/ncbienv.hpp>
#include <corelib/ncbiargs.hpp>

#include <objmgr/object_manager.hpp>
#include <objmgr/util/sequence.hpp>

#include <objects/seqalign/Seq_align.hpp>
#include <objects/seqalign/Seq_align_set.hpp>

#include <objtools/blast/seqdb_reader/seqdb.hpp>
#include <objtools/blast/seqdb_reader/seqdbcommon.hpp>

#include <algo/blast/api/sseqloc.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/uniform_search.hpp>
#include <algo/blast/api/blast_types.hpp>
#include <algo/blast/api/blast_aux.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/api/blast_nucl_options.hpp>
#include <algo/blast/api/blast_prot_options.hpp>
#include <algo/blast/api/remote_blast.hpp>
#include "blast_app_util.hpp"

#include <algo/blast/blastinput/blast_input.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>
#include <algo/blast/blastinput/blast_args.hpp>

#include <algorithm>
#include <string>
#include <cmath>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <forward_list>
#include <chrono>

USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(objects);


/// Corresponds to ranks in NCBI taxonomy database.
enum ETaxRank {
  eRank_No_Rank = -1,

  eRank_Subspecies = 0,
  eRank_Forma = eRank_Subspecies,
  eRank_Varietas = eRank_Subspecies,
  eRank_Species,
  eRank_Species_Subgroup,
  eRank_Species_Group,
  eRank_Subgenus,
  eRank_Genus,
  eRank_Subtribe,
  eRank_Tribe,
  eRank_Subfamily,
  eRank_Family,
  eRank_Superfamily,
  eRank_Parvorder,
  eRank_Infraorder,
  eRank_Suborder,
  eRank_Order,
  eRank_Superorder,
  eRank_Infraclass,
  eRank_Subclass,
  eRank_Class,
  eRank_Superclass,
  eRank_Subphylum,
  eRank_Phylum,
  eRank_Superphylum,
  eRank_Subkingdom,
  eRank_Kingdom,
  eRank_Superkingdom,

  eRank_All // For filter purposes
};


/// Converts a string-formatted taxon rank to ETaxRank value.
/// @param name Taxon rank name [in]
ETaxRank RankNameToEnum(const string &name) {
  if      (name == "species")          return eRank_Species;
  else if (name == "subspecies")       return eRank_Subspecies;
  else if (name == "forma")            return eRank_Forma;
  else if (name == "varietas")         return eRank_Varietas;
  else if (name == "species subgroup") return eRank_Species_Subgroup;
  else if (name == "species group")    return eRank_Species_Group;
  else if (name == "no rank")          return eRank_No_Rank;
  else if (name == "subgenus")         return eRank_Subgenus;
  else if (name == "genus")            return eRank_Genus;
  else if (name == "subtribe")         return eRank_Subtribe;
  else if (name == "tribe")            return eRank_Tribe;
  else if (name == "subfamily")        return eRank_Subfamily;
  else if (name == "family")           return eRank_Family;
  else if (name == "superfamily")      return eRank_Superfamily;
  else if (name == "parvorder")        return eRank_Parvorder;
  else if (name == "infraorder")       return eRank_Infraorder;
  else if (name == "suborder")         return eRank_Suborder;
  else if (name == "order")            return eRank_Order;
  else if (name == "superorder")       return eRank_Superorder;
  else if (name == "infraclass")       return eRank_Infraclass;
  else if (name == "subclass")         return eRank_Subclass;
  else if (name == "class")            return eRank_Class;
  else if (name == "superclass")       return eRank_Superclass;
  else if (name == "subphylum")        return eRank_Subphylum;
  else if (name == "phylum")           return eRank_Phylum;
  else if (name == "superphylum")      return eRank_Superphylum;
  else if (name == "subkingdom")       return eRank_Subkingdom;
  else if (name == "kingdom")          return eRank_Kingdom;
  else if (name == "superkingdom")     return eRank_Superkingdom;
  else if (name == "all")              return eRank_All;
  cerr << "Unknown tax rank " << name << endl;
  return eRank_No_Rank; // Default
}


/// Tokenize a std::string and place the output in a premade string container. Must
/// be a STL or STL-compatible container type that implements insert(position, value).
/// @param in String to tokenize [in]
/// @param delimiters String containing characters to use as token delimiter(s) [in]
/// @param out Already created string container to place results in [out]
template <class containerT>
void TokenizeString(const string &in, const string &delimiters, containerT &out) {
  using valueT = typename containerT::value_type;
  using sizeT = typename containerT::size_type;

  string::size_type curPosition, lastPosition = 0;

  while (lastPosition < in.length()) {
    curPosition = in.find_first_of(delimiters, lastPosition);
    if (curPosition == string::npos)
      curPosition = in.length();

    if (curPosition != lastPosition)
      out.insert(out.end(),
        valueT(in.data() + lastPosition, (sizeT)(curPosition - lastPosition)));

    lastPosition = curPosition + 1;
  }
}


/// Tokenize a std::string and place the output in a premade integer container. Must
/// be a STL or STL-compatible container type that implements insert(position, value).
/// @param in String to tokenize [in]
/// @param delimiters String containing characters to use as token delimiter(s) [in]
/// @param out Already created integer container to place results in [out]
template <class containerT>
void TokenizeStrToInt(const string &in, const string &delimiters, containerT &out) {
  vector<string> tokens;
  TokenizeString(in, delimiters, tokens);

  for (string &substr: tokens) {
    int result;
    try {
      result = stoi(substr);
    }
    catch (...) {
      continue; // Could not convert substring to int, skip
    }

    out.insert(out.end(), result);
  }
}


/////////////////////////////////////////////////////////////////////////////
//  CAlienG2::


class CAlienG2 : public CNcbiApplication {
public:
  typedef unordered_set<int> TGroup;
  typedef tuple<CRef<CSeq_align>, CRef<CSeq_align>, /* Group 1, group 2 */
                CSeq_id_Handle, CSeq_id_Handle,
                CBioseq_Handle, CBioseq_Handle,
                double, double, double, /* Alien index, HGT index, Ratio */
                int, int /* Tax ID 1, Tax ID 2 */> THgtData;

private:
  virtual void Init();
  virtual int Run();
  virtual void Exit();

  void SetupFiltersAndGroups();
  void ProcessBlastOptions(CRef<CBlastOptionsHandle> optsHandle);
  bool ProcessTaxonomyDump();

  bool IsTaxonInGroup(int taxId, TGroup &group);
  int GetCommonAncestor(int taxId1, int taxId2, ETaxRank *highestRank = nullptr);

  void EvaluateResults(CRef<CSearchResultSet> results, CRef<CScope> scope);
  void CreateReport(CNcbiOstream &out);

  CRef<CBlastQueryVector> x_ExtractQueries(bool isProtein);
  CRef<CBlastSearchQuery> x_BuildQueryFromPssm(const CPssmWithParameters &pssm);
  SSeqLoc x_QueryBioseqToSSeqLoc(const CBioseq& bioseq, CRef<CScope> scope);

  vector<int> parentTaxId;
  vector<ETaxRank> taxRank;
  vector<string> sciName;
  unordered_multimap<string, int> taxNameId;

  unordered_set<int> filteredGIs;
  unordered_set<string> filteredAccs;
  TGroup filteredTaxa,
         group1Taxa,
         group2Taxa;
  double bitscoreThreshold,
         eValueThreshold,
         coverageThreshold,
         alienIndexMin,
         alienIndexMax,
         hgtIndexMin,
         hgtIndexMax,
         scoreRatioMin,
         scoreRatioMax;
  ETaxRank minAncestorRank;
  bool group1First,
       group2First;

  list<THgtData> hgtData;

  CRef<CLocalBlast> lclBlast;
  CRef<CRemoteBlast> rmtBlast;
  CRef<CBlastScopeSource> queryScopeSource;
};


/////////////////////////////////////////////////////////////////////////////
// Initialize parameters/arguments
//
// Command-line arguments are checked first. If they are not set, then config
// parameters are checked.


// Define config parameters
NCBI_PARAM_DECL(string, AlienG2, task);
NCBI_PARAM_DEF( string, AlienG2, task, "report");
NCBI_PARAM_DECL(int,    AlienG2, num_threads);
NCBI_PARAM_DEF( int,    AlienG2, num_threads, (int)CThreadable::kMinNumThreads);
NCBI_PARAM_DECL(string, AlienG2, database);
NCBI_PARAM_DEF( string, AlienG2, database, "nr");
NCBI_PARAM_DECL(string, AlienG2, taxdump);
NCBI_PARAM_DEF( string, AlienG2, taxdump, "");
NCBI_PARAM_DECL(string, AlienG2, in);
NCBI_PARAM_DEF( string, AlienG2, in, "stdin");
NCBI_PARAM_DECL(string, AlienG2, out);
NCBI_PARAM_DEF( string, AlienG2, out, "stdout");
typedef NCBI_PARAM_TYPE(AlienG2, task)        TParam_Task;
typedef NCBI_PARAM_TYPE(AlienG2, num_threads) TParam_NumThreads;
typedef NCBI_PARAM_TYPE(AlienG2, database)    TParam_Database;
typedef NCBI_PARAM_TYPE(AlienG2, taxdump)     TParam_TaxDump;
typedef NCBI_PARAM_TYPE(AlienG2, in)          TParam_In;
typedef NCBI_PARAM_TYPE(AlienG2, out)         TParam_Out;

NCBI_PARAM_DECL(int,    AlienG2, max_target_seqs);
NCBI_PARAM_DEF( int,    AlienG2, max_target_seqs, 0);
NCBI_PARAM_DECL(int,    AlienG2, penalty);
NCBI_PARAM_DEF( int,    AlienG2, penalty, 0);
NCBI_PARAM_DECL(int,    AlienG2, reward);
NCBI_PARAM_DEF( int,    AlienG2, reward, 0);
NCBI_PARAM_DECL(string, AlienG2, matrix);
NCBI_PARAM_DEF( string, AlienG2, matrix, "");
NCBI_PARAM_DECL(string, AlienG2, evalue);
NCBI_PARAM_DEF( string, AlienG2, evalue, "0");
NCBI_PARAM_DECL(string, AlienG2, bitscore);
NCBI_PARAM_DEF( string, AlienG2, bitscore, "0");
NCBI_PARAM_DECL(string, AlienG2, coverage);
NCBI_PARAM_DEF( string, AlienG2, coverage, "0");
typedef NCBI_PARAM_TYPE(AlienG2, max_target_seqs) TParam_MaxTargetSeqs;
typedef NCBI_PARAM_TYPE(AlienG2, penalty)         TParam_Penalty;
typedef NCBI_PARAM_TYPE(AlienG2, reward)          TParam_Reward;
typedef NCBI_PARAM_TYPE(AlienG2, matrix)          TParam_Matrix;
typedef NCBI_PARAM_TYPE(AlienG2, evalue)          TParam_EValue;
typedef NCBI_PARAM_TYPE(AlienG2, bitscore)        TParam_Bitscore;
typedef NCBI_PARAM_TYPE(AlienG2, coverage)        TParam_Coverage;

NCBI_PARAM_DECL(string, AlienG2, filter);
NCBI_PARAM_DEF( string, AlienG2, filter, "");
NCBI_PARAM_DECL(string, AlienG2, filter_ids);
NCBI_PARAM_DEF( string, AlienG2, filter_ids, "");
NCBI_PARAM_DECL(string, AlienG2, filter_accs);
NCBI_PARAM_DEF( string, AlienG2, filter_accs, "");
NCBI_PARAM_DECL(string, AlienG2, filter_gis);
NCBI_PARAM_DEF( string, AlienG2, filter_gis, "");
NCBI_PARAM_DECL(string, AlienG2, group1);
NCBI_PARAM_DEF( string, AlienG2, group1, "");
NCBI_PARAM_DECL(string, AlienG2, group1_ids);
NCBI_PARAM_DEF( string, AlienG2, group1_ids, "");
NCBI_PARAM_DECL(string, AlienG2, group2);
NCBI_PARAM_DEF( string, AlienG2, group2, "");
NCBI_PARAM_DECL(string, AlienG2, group2_ids);
NCBI_PARAM_DEF( string, AlienG2, group2_ids, "");
typedef NCBI_PARAM_TYPE(AlienG2, filter)      TParam_FilterNames;
typedef NCBI_PARAM_TYPE(AlienG2, filter_ids)  TParam_FilterIDs;
typedef NCBI_PARAM_TYPE(AlienG2, filter_accs) TParam_FilterAccs;
typedef NCBI_PARAM_TYPE(AlienG2, filter_gis)  TParam_FilterGIs;
typedef NCBI_PARAM_TYPE(AlienG2, group1)      TParam_Group1Names;
typedef NCBI_PARAM_TYPE(AlienG2, group1_ids)  TParam_Group1IDs;
typedef NCBI_PARAM_TYPE(AlienG2, group2)      TParam_Group2Names;
typedef NCBI_PARAM_TYPE(AlienG2, group2_ids)  TParam_Group2IDs;

NCBI_PARAM_DECL(string, AlienG2, alien_min);
NCBI_PARAM_DEF( string, AlienG2, alien_min, "0");
NCBI_PARAM_DECL(string, AlienG2, alien_max);
NCBI_PARAM_DEF( string, AlienG2, alien_max, "0");
NCBI_PARAM_DECL(string, AlienG2, hgt_min);
NCBI_PARAM_DEF( string, AlienG2, hgt_min, "0");
NCBI_PARAM_DECL(string, AlienG2, hgt_max);
NCBI_PARAM_DEF( string, AlienG2, hgt_max, "0");
NCBI_PARAM_DECL(string, AlienG2, ratio_min);
NCBI_PARAM_DEF( string, AlienG2, ratio_min, "0");
NCBI_PARAM_DECL(string, AlienG2, ratio_max);
NCBI_PARAM_DEF( string, AlienG2, ratio_max, "0");
NCBI_PARAM_DECL(string, AlienG2, ancestor_rank);
NCBI_PARAM_DEF( string, AlienG2, ancestor_rank, "no rank");
NCBI_PARAM_DECL(string, AlienG2, search_order);
NCBI_PARAM_DEF( string, AlienG2, search_order, "any");
typedef NCBI_PARAM_TYPE(AlienG2, alien_min)     TParam_AlienIndexMin;
typedef NCBI_PARAM_TYPE(AlienG2, alien_max)     TParam_AlienIndexMax;
typedef NCBI_PARAM_TYPE(AlienG2, hgt_min)       TParam_HGTIndexMin;
typedef NCBI_PARAM_TYPE(AlienG2, hgt_max)       TParam_HGTIndexMax;
typedef NCBI_PARAM_TYPE(AlienG2, ratio_min)     TParam_ScoreRatioMin;
typedef NCBI_PARAM_TYPE(AlienG2, ratio_max)     TParam_ScoreRatioMax;
typedef NCBI_PARAM_TYPE(AlienG2, ancestor_rank) TParam_AncestorRank;
typedef NCBI_PARAM_TYPE(AlienG2, search_order)  TParam_SearchOrder;


/// Initialize command line arguments. Default values taken from config parameters.
void CAlienG2::Init() {
  auto_ptr<CArgDescriptions> argDesc(new CArgDescriptions);

  // Specify USAGE context
  argDesc->SetUsageContext(GetArguments().GetProgramBasename(),
    "AlienG2 HGT Analyzer");

  // ** TODO Need to test all programs to ensure they work
  TParam_Task pTask;
  argDesc->AddDefaultKey(
    "task", "task_name",
    "One of blastn, megablast, dc-megablast, blastp, blastx, tblastn, "
    "tblastx, rpsblast, rpstblastn, psiblast, psitblastn, deltablast"
    /*", phiblastp, phiblastn, vecscreen"*/
    ", report (default)",
    CArgDescriptions::eString, pTask.Get());
  argDesc->SetConstraint(
    "task", &(*new CArgAllow_Strings,
    "blastn", "megablast", "dc-megablast", "blastp", "blastx", "tblastn",
    "tblastx", "rpsblast", "rpstblastn", "psiblast", "psitblastn",
    "deltablast", "report", ""/*, "phiblastp", "phiblastn", "vecscreen"*/));

  TParam_NumThreads pNumThreads;
  argDesc->AddDefaultKey
    ("num_threads", "int_value", "Number of BLAST search threads to spawn",
    CArgDescriptions::eInteger, to_string(pNumThreads.Get()));

  TParam_Database pDatabase;
  argDesc->AddDefaultKey
    ("db", "database_name", "Name of database (such as nr, nt) for search",
    CArgDescriptions::eString, pDatabase.Get());

  TParam_TaxDump pTaxDump;
  argDesc->AddDefaultKey
    ("taxdump", "directory", "Directory containing NCBI taxdump files",
    CArgDescriptions::eString, pTaxDump.Get());

  TParam_In pInput;
  argDesc->AddDefaultKey("in", "input_file",
    "A report file in BLAST Archive format (outfmt 11), or a FASTA file with "
    "the query with -task= BLAST",
    CArgDescriptions::eInputFile, pInput.Get());

  TParam_Out pOutput;
  argDesc->AddDefaultKey("out", "output_file",
    "The output file", CArgDescriptions::eOutputFile, pOutput.Get());

  TParam_MaxTargetSeqs pMaxTargetSeqs;
  argDesc->AddDefaultKey("max_target_seqs", "num_sequences",
    "Maximum number of aligned sequences to keep",
    CArgDescriptions::eInteger, to_string(pMaxTargetSeqs.Get()));

  TParam_Penalty pPenalty;
  argDesc->AddDefaultKey("penalty", "nucl_mismatch_penalty",
    "Penalty score for a mismatch (with -task= nucleotide BLAST only)",
    CArgDescriptions::eInteger, to_string(pPenalty.Get()));

  TParam_Reward pReward;
  argDesc->AddDefaultKey("reward", "nucl_match_reward",
    "Reward score for a match (with -task= nucleotide BLAST only)",
    CArgDescriptions::eInteger, to_string(pReward.Get()));

  TParam_Matrix pMatrix;
  argDesc->AddDefaultKey("matrix", "prot_matrix_name",
    "Scoring matrix name (with -task= protein BLAST only)",
    CArgDescriptions::eString, pMatrix.Get());

  TParam_EValue pEValue;
  argDesc->AddDefaultKey("evalue", "evalue_max",
    "E-value maximum threshold for saving hits (0 = none)",
    CArgDescriptions::eDouble, pEValue.Get());

  TParam_Bitscore pBitscore;
  argDesc->AddDefaultKey("bitscore", "bitscore_min",
    "Bitscore minimum threshold for saving hits (0 = none)",
    CArgDescriptions::eDouble, pBitscore.Get());

  TParam_Coverage pCoverage;
  argDesc->AddDefaultKey("coverage", "coverage_min",
    "Percent coverage (align length ratio) threshold for saving hits",
    CArgDescriptions::eDouble, pCoverage.Get());

  TParam_FilterNames pFilterNames;
  argDesc->AddDefaultKey("filter", "taxon_names",
    "Taxon names (in taxdump) to filter from results, separated by semicolons "
    "(no trailing space)\n"
    "NOTE: A name can match multiple taxa (ex. \"environmental samples\")",
    CArgDescriptions::eString, pFilterNames.Get());

  TParam_FilterIDs pFilterIDs;
  argDesc->AddDefaultKey("filter_ids", "taxon_id_list",
    "Taxon IDs to filter from results, separated by semicolons",
    CArgDescriptions::eString, pFilterIDs.Get());

  TParam_FilterAccs pFilterAccs;
  argDesc->AddDefaultKey("filter_accs", "accession_list",
    "Sequence accession numbers to filter from results, separated by semicolons "
    "(no trailing space)",
    CArgDescriptions::eString, pFilterAccs.Get());

  TParam_FilterGIs pFilterGIs;
  argDesc->AddDefaultKey("filter_gis", "gi_list",
    "Sequence GIs to filter from results, separated by semicolons",
    CArgDescriptions::eString, pFilterGIs.Get());

  TParam_Group1Names pGroup1Names;
  argDesc->AddDefaultKey("group1", "taxon_name_list",
    "Taxon names (in taxdump) part of HGT group 1, separated by semicolons "
    "(no trailing space)\n"
    "NOTE: A name can match multiple taxa (ex. \"environmental samples\")",
    CArgDescriptions::eString, pGroup1Names.Get());

  TParam_Group1IDs pGroup1IDs;
  argDesc->AddDefaultKey("group1_ids", "taxon_id_list",
    "Taxon IDs part of HGT group 1, separated by semicolons",
    CArgDescriptions::eString, pGroup1IDs.Get());

  TParam_Group2Names pGroup2Names;
  argDesc->AddDefaultKey("group2", "taxon_name_list",
    "Taxon names (in taxdump) part of HGT group 2, separated by semicolons "
    "(no trailing space)\n"
    "NOTE: A name can match multiple taxa (ex. \"environmental samples\")",
    CArgDescriptions::eString, pGroup2Names.Get());

  TParam_Group2IDs pGroup2IDs;
  argDesc->AddDefaultKey("group2_ids", "taxon_id_list",
    "Taxon IDs part of HGT group 2, separated by semicolons",
    CArgDescriptions::eString, pGroup2IDs.Get());

  TParam_AlienIndexMin pAlienIndexMin;
  argDesc->AddDefaultKey("alien_min", "score_min",
    "Alien index minimum threshold for saving hit pairs (0 = none)",
    CArgDescriptions::eDouble, pAlienIndexMin.Get());

  TParam_AlienIndexMin pAlienIndexMax;
  argDesc->AddDefaultKey("alien_max", "score_max",
    "Alien index maximum threshold for saving hit pairs (0 = none)",
    CArgDescriptions::eDouble, pAlienIndexMax.Get());

  TParam_HGTIndexMin pHGTIndexMin;
  argDesc->AddDefaultKey("hgt_min", "score_min",
    "HGT index minimum threshold for saving hit pairs (0 = none)",
    CArgDescriptions::eDouble, pHGTIndexMin.Get());

  TParam_HGTIndexMin pHGTIndexMax;
  argDesc->AddDefaultKey("hgt_max", "score_max",
    "HGT index maximum threshold for saving hit pairs (0 = none)",
    CArgDescriptions::eDouble, pHGTIndexMax.Get());

  TParam_ScoreRatioMin pScoreRatioMin;
  argDesc->AddDefaultKey("ratio_min", "score_min",
    "Group 1:2 bitscore ratio minimum threshold for saving hit pairs (0 = none)",
    CArgDescriptions::eDouble, pScoreRatioMin.Get());

  TParam_ScoreRatioMax pScoreRatioMax;
  argDesc->AddDefaultKey("ratio_max", "score_max",
    "Group 1:2 bitscore ratio maximum threshold for saving hit pairs (0 = none)",
    CArgDescriptions::eDouble, pScoreRatioMax.Get());

  TParam_AncestorRank pAncestorRank;
  argDesc->AddDefaultKey("ancestor_rank", "taxon_rank",
    "Filter results whose nearest common taxon with the better-scoring group is "
    "not at least this rank (no rank = no threshold, all = no common superkingdom)",
    CArgDescriptions::eString, pAncestorRank.Get());
  argDesc->SetConstraint("ancestor_rank", &(*new CArgAllow_Strings,
    "no rank", "subspecies", "forma", "varietas", "species", "species subgroup",
    "species group", "subgenus", "genus", "subtribe", "tribe", "subfamily",
    "family", "superfamily", "parvorder", "infraorder", "suborder", "order",
    "superorder", "infraclass", "subclass", "class", "superclass", "subphylum",
    "phylum", "superphylum", "subkingdom", "kingdom", "superkingdom", "all"));

  TParam_SearchOrder pSearchOrder;
  argDesc->AddDefaultKey("search_order", "group_priority",
    "\"any\" = Either group can come first, \"group1\" = Group 1 must come "
    "first, \"group2\" = Group 2 must come first",
    CArgDescriptions::eString, pSearchOrder.Get());
  argDesc->SetConstraint("search_order",
    &(*new CArgAllow_Strings, "any", "group1", "group2"));

  // Setup arg.descriptions for this application
  SetupArgDescriptions(argDesc.release());
}


/// Parse NCBI taxonomy dump files to get taxonomic information.
bool CAlienG2::ProcessTaxonomyDump() {
  const CArgs &args = GetArgs();

  string taxDumpDir = args["taxdump"].AsString();
  if (!taxDumpDir.empty() && taxDumpDir.back() != '/' && taxDumpDir.back() != '\\')
    taxDumpDir += taxDumpDir.find('\\') != string::npos ? '\\' : '/';

  // Format is text, where rows are separated by \n, fields separated by \t|\t

  // Parse nodes.dmp and populate maps
  CNcbiIfstream taxNodes(taxDumpDir + "nodes.dmp");
  if (taxNodes.is_open()) {
    string l;
    char *strtolEnd;

    unordered_set<string> filtered;
    TokenizeString(args["filter"].AsString(), ";", filtered);
    bool loadNames = !(filtered.empty() && args["group1"].AsString().empty() &&
                       args["group2"].AsString().empty()),
         filterUncultured = filtered.count("environmental samples") != 0;
    filtered.clear();

    // Reserve space in vectors to prevent constantly resizing
    const int maxNodes = 2250000; // Set to >= # of lines in nodes+merged * 1.3
    parentTaxId.resize(maxNodes);
    taxRank.resize(maxNodes);
    sciName.resize(maxNodes);
    // Reserve space in maps to prevent constantly rehashing
    if (loadNames)
      taxNameId.reserve(3750000); // Set to >= # of names + uniques * 1.3
    unordered_map<int, int> merged;
    merged.reserve(65000); // Set to >= number of lines in merged * 1.3

    while (!taxNodes.eof()) {
      getline(taxNodes, l, '\n');
      if (l.empty())
        continue;

      // Get first 3 fields: tax_id, parent tax_id, rank
      unsigned int taxId       = strtol(l.data(), &strtolEnd, 10),
                   pTaxIdStart = l.find_first_of('\t', strtolEnd - l.data()) + 3;

      if (taxId >= parentTaxId.size()) {
        parentTaxId.resize(taxId + 1000);
        taxRank.resize(taxId + 1000);
        sciName.resize(taxId + 1000);
      }

      parentTaxId[taxId] = strtol(l.data() + pTaxIdStart, &strtolEnd, 10);

      int rankStart = l.find_first_of('\t', strtolEnd - l.data()) + 3,
          rankEnd   = l.find_first_of('\t', rankStart);
      taxRank[taxId] = RankNameToEnum(l.substr(rankStart, rankEnd - rankStart));

      // If we are filtering taxa named "environmental samples", and this
      // taxon's comments say "uncultured", then add it to ID filter
      if (filterUncultured && filteredTaxa.count(taxId) == 0) {
        int cEnd   = l.find_last_of('\t'),
            cStart = l.find_last_of('\t', cEnd - 1) + 1;
        if (cEnd - cStart >= 2 && strncmp(l.data() + cStart, "un",
            2) == 0)
          filteredTaxa.insert(taxId);
      }
    }
    taxNodes.close();

    // Parse merged.dmp and add entries from merged nodes to maps for compatibility
    CNcbiIfstream mergedNodes(taxDumpDir + "merged.dmp");
    if (mergedNodes.is_open()) {
      while (!mergedNodes.eof()) {
        getline(mergedNodes, l, '\n');
        if (l.empty())
          continue;

        // Get first 2 fields: old_tax_id, new_tax_id
        unsigned int oTaxId = strtol(l.data(), &strtolEnd, 10),
                     nTaxId = strtol(l.data() +
                       l.find_first_of('\t', strtolEnd - l.data()) + 3, nullptr, 10);

        if (oTaxId >= parentTaxId.size()) {
          parentTaxId.resize(oTaxId + 1000);
          taxRank.resize(oTaxId + 1000);
          sciName.resize(oTaxId + 1000);
        }

        parentTaxId[oTaxId] = parentTaxId[nTaxId];
        taxRank[oTaxId] = taxRank[nTaxId];

        // Cache this so we can copy names later
        merged[nTaxId] = oTaxId;
      }
      mergedNodes.close();
    }

    // Parse names.dmp and populate names map
    CNcbiIfstream nodeNames(taxDumpDir + "names.dmp");
    if (nodeNames.is_open()) {
      while (!nodeNames.eof()) {
        getline(nodeNames, l, '\n');
        if (l.empty())
          continue;

        // Get first 3 fields: tax_id, name, unique_name
        unsigned int taxId = strtol(l.data(), &strtolEnd, 10);
        int nameStart = l.find_first_of('\t', strtolEnd - l.data()) + 3,
            nameEnd = l.find_first_of('\t', nameStart),
            uNameStart = nameEnd + 3,
            uNameEnd = l.find_first_of('\t', uNameStart);

        if (taxId >= sciName.size())
          sciName.resize(taxId + 1000);

        // If type is scientific name, cache the name (used in making report)
        if (nameStart != nameEnd && sciName[taxId].empty() &&
            strncmp(l.data() + uNameEnd + 3, "sc", 2) == 0) {
          sciName[taxId].assign(l.data() + nameStart, nameEnd - nameStart);

          if (merged.count(taxId) != 0)
            sciName[merged[taxId]].assign(l.data() + nameStart, nameEnd - nameStart);
        }

        if (loadNames) {
          // Change to lowercase for case-insensitive matching
          transform(l.begin() + nameStart, l.begin() + uNameEnd,
            l.begin() + nameStart, ::tolower);

          // Emplace name (and unique name) into names map
          if (nameStart != nameEnd)
            taxNameId.emplace(piecewise_construct,
              forward_as_tuple(l.data() + nameStart, nameEnd - nameStart),
              forward_as_tuple(taxId));
          if (uNameStart != uNameEnd)
            taxNameId.emplace(piecewise_construct,
              forward_as_tuple(l.data() + uNameStart, uNameEnd - uNameStart),
              forward_as_tuple(taxId));

          if (merged.count(taxId) != 0) {
            // Copy names to original merged node
            if (nameStart != nameEnd)
              taxNameId.emplace(piecewise_construct,
                forward_as_tuple(l.data() + nameStart, nameEnd - nameStart),
                forward_as_tuple(merged[taxId]));
            if (uNameStart != uNameEnd)
              taxNameId.emplace(piecewise_construct,
                forward_as_tuple(l.data() + uNameStart, uNameEnd - uNameStart),
                forward_as_tuple(merged[taxId]));
          }
        }
      }
      nodeNames.close();
    }
    else
      return false;

    return true;
  }
  else
    return false;
}


/// Setup search filters and HGT groups based upon arguments/parameters.
void CAlienG2::SetupFiltersAndGroups() {
  const CArgs &args = GetArgs();

  vector<string> names;

  // Populate sequence filter
  string sFilterAccs = args["filter_accs"].AsString(),
         sFilterGIs  = args["filter_gis"].AsString();
  TokenizeString(sFilterAccs, ";", filteredAccs);
  TokenizeStrToInt(sFilterGIs, ";", filteredGIs);

  // Populate taxon filter
  string sFilterNames = args["filter"].AsString(),
         sFilterIds   = args["filter_ids"].AsString();
  transform(sFilterNames.begin(), sFilterNames.end(), sFilterNames.begin(),
    ::tolower); // Change to lowercase for case-insensitive matching
  TokenizeString(sFilterNames, ";", names);
  for (auto &name: names) {
    auto range = taxNameId.equal_range(name);
    for (auto curName = range.first; curName != range.second; ++curName)
      filteredTaxa.insert(curName->second);
  }
  names.clear();
  TokenizeStrToInt(sFilterIds, ";", filteredTaxa);

  // Populate HGT group 1 criteria
  string sGroup1Names = args["group1"].AsString(),
         sGroup1Ids   = args["group1_ids"].AsString();
  transform(sGroup1Names.begin(), sGroup1Names.end(), sGroup1Names.begin(),
    ::tolower);
  TokenizeString(sGroup1Names, ";", names);
  for (auto &name: names) {
    auto range = taxNameId.equal_range(name);
    for (auto curName = range.first; curName != range.second; ++curName)
      group1Taxa.insert(curName->second);
  }
  names.clear();
  TokenizeStrToInt(sGroup1Ids, ";", group1Taxa);

  // Populate HGT group 2 criteria
  string sGroup2Names = args["group2"].AsString(),
         sGroup2Ids   = args["group2_ids"].AsString();
  if (sGroup2Names == "*")
    group2Taxa.insert(1);
  else {
    transform(sGroup2Names.begin(), sGroup2Names.end(), sGroup2Names.begin(),
      ::tolower);
    TokenizeString(sGroup2Names, ";", names);
    for (auto &name: names) {
      auto range = taxNameId.equal_range(name);
      for (auto curName = range.first; curName != range.second; ++curName)
        group2Taxa.insert(curName->second);
    }
    names.clear();
    TokenizeStrToInt(sGroup2Ids, ";", group2Taxa);
  }

  // ** TODO Remove this when null group 2 is supported
  // HGT Group 1 cannot contain root node
  if (group1Taxa.count(1) != 0)
    group1Taxa.erase(1);
  // If empty, add root node to HGT group 2 criteria
  if (group2Taxa.empty())
    group2Taxa.insert(1);

  // Unload full tax name map, we are done using it
  taxNameId.clear();
}

/////////////////////////////////////////////////////////////////////////////
//  Query initialization functions taken from blast_formatter.cpp (NCBI)


/// Package a scope and Seq-loc into a SSeqLoc from a Bioseq
/// @param bioseq Bioseq to inspect [in]
/// @param scope Scope object to add the sequence data to [in|out]
SSeqLoc CAlienG2::x_QueryBioseqToSSeqLoc(const CBioseq& bioseq, CRef<CScope> scope) {
  static bool firstTime = true;
  _ASSERT(scope);

  if (!HasRawSequenceData(bioseq) && firstTime) {
    _ASSERT(queryScopeSource);
    queryScopeSource->AddDataLoaders(scope);
    firstTime = false;
  }
  else
    scope->AddBioseq(bioseq);

  CRef<CSeq_loc> seqloc(new CSeq_loc);
  seqloc->SetWhole().Assign(*bioseq.GetFirstId());
  return SSeqLoc(seqloc, scope);
}


/// Build the query from a PSSM
/// @param pssm PSSM to inspect [in]
CRef<CBlastSearchQuery>
CAlienG2::x_BuildQueryFromPssm(const CPssmWithParameters &pssm) {
  if (!pssm.HasQuery())
    throw runtime_error("PSSM has no query");

  CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
  const CSeq_entry& seqEntry = pssm.GetQuery();
  if (!seqEntry.IsSeq())
    throw runtime_error("Cannot have multiple queries in a PSSM");

  SSeqLoc ssl = x_QueryBioseqToSSeqLoc(seqEntry.GetSeq(), scope);
  CRef<CBlastSearchQuery> retval;
  retval.Reset(new CBlastSearchQuery(*ssl.seqloc, *ssl.scope));
  _ASSERT(ssl.scope.GetPointer() == scope.GetPointer());
  return retval;
}


/// Extracts queries from report file and adds them to scope
/// @param isProtein Are the queries protein sequences? [in]
CRef<CBlastQueryVector> CAlienG2::x_ExtractQueries(bool isProtein) {
  CRef<CBlast4_queries> b4_queries = rmtBlast->GetQueries();
  _ASSERT(b4_queries);
  const size_t kNumQueries = b4_queries->GetNumQueries();

  CRef<CBlastQueryVector> retval(new CBlastQueryVector);

  SDataLoaderConfig dlconfig(isProtein, SDataLoaderConfig::eUseNoDataLoaders);
  dlconfig.OptimizeForWholeLargeSequenceRetrieval(false);
  queryScopeSource.Reset(new CBlastScopeSource(dlconfig));

  if (b4_queries->IsPssm())
    retval->AddQuery(x_BuildQueryFromPssm(b4_queries->GetPssm()));
  else if (b4_queries->IsSeq_loc_list()) {
    CRef<CScope> scope = queryScopeSource->NewScope();
    ITERATE(CBlast4_queries::TSeq_loc_list, seqloc,
        b4_queries->GetSeq_loc_list()) {
      _ASSERT(!(*seqloc)->GetId()->IsLocal());
      CRef<CBlastSearchQuery> query(new CBlastSearchQuery(**seqloc,
        *scope));
      retval->AddQuery(query);
    }
  }
  else if (b4_queries->IsBioseq_set()) {
    CTypeConstIterator<CBioseq> itr(ConstBegin(b4_queries->GetBioseq_set(),
      eDetectLoops));
    CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
    for (; itr; ++itr) {
      SSeqLoc ssl = x_QueryBioseqToSSeqLoc(*itr, scope);
      CRef<CBlastSearchQuery> query(new CBlastSearchQuery(*ssl.seqloc,
        *ssl.scope));
      retval->AddQuery(query);
    }
  }

  (void)kNumQueries; // Eliminate compiler warning
  _ASSERT(kNumQueries == retval->size());
  return retval;
}


//
/////////////////////////////////////////////////////////////////////////////


/// Modify BLAST options from defaults based upon command-line args or config.
/// @param optsHandle already created CBlastOptionsHandle to modify [in]
void CAlienG2::ProcessBlastOptions(CRef<CBlastOptionsHandle> optsHandle) {
  const CArgs &args = GetArgs();

  // Expect value is a supported option for all flavors of BLAST
  if (args["evalue"].AsDouble())
    optsHandle->SetEvalueThreshold(args["evalue"].AsDouble());

  if (args["max_target_seqs"].AsInteger() > 0)
    optsHandle->SetHitlistSize(args["max_target_seqs"].AsInteger());

  // The first branch is used if the program is blastn or a flavor of megablast
  // as reward and penalty is a valid option
  // The second branch is used for all other programs except rpsblast as matrix
  // is a valid option for blastp and other programdeltablasts that perform
  // protein-protein comparisons
  if (CBlastNucleotideOptionsHandle* nuclHandle =
      dynamic_cast<CBlastNucleotideOptionsHandle*>(&*optsHandle)) {
    if (args["reward"].AsInteger())
      nuclHandle->SetMatchReward(args["reward"].AsInteger());
    if (args["penalty"].AsInteger())
      nuclHandle->SetMismatchPenalty(args["penalty"].AsInteger());
  }
  else if (CBlastProteinOptionsHandle* protHandle =
      dynamic_cast<CBlastProteinOptionsHandle*>(&*optsHandle)) {
    if (args["matrix"] && !args["matrix"].AsString().empty())
      protHandle->SetMatrixName(args["matrix"].AsString().c_str());
  }

  return;
}


/// Returns true if a taxon or one of its parents matches criteria specified
/// in @group.
/// @param taxId Taxon ID to query [in]
/// @param group Taxon group to test against [in/out]
bool CAlienG2::IsTaxonInGroup(int taxId, TGroup &group) {
  if (taxId < 1 || group.empty()) // Invalid ID or no criteria
    return false;
  if (group.count(1) != 0)
    return true; // Group containing root node contains all taxa

  // Test if this taxId or any of its ancestors meet criteria
  for (int i = taxId; i > 1; i = parentTaxId[i]) {
    if (group.count(i) != 0) {
      if (i != taxId)
        group.insert(taxId);
      return true;
    }
  }

  return false;
}


/// Get nearest common ancestor for two taxa.
/// @param taxId1 First taxon ID to query [in]
/// @param taxId2 Second taxon ID to query [in]
/// @param highestRank Highest taxon ranking of ancestor [out] (optional)
int CAlienG2::GetCommonAncestor(int taxId1, int taxId2, ETaxRank *highestRank) {
  if (taxId1 == taxId2) {
    if (highestRank != nullptr)
      *highestRank = taxRank[taxId1];
    return taxId1;
  }
  if (taxId1 <= 1 || taxId2 <= 1) { // Invalid taxId
    if (highestRank != nullptr)
      *highestRank = taxRank[1]; // eRank_No_Rank
    return 1;
  }

  // Many intermediary taxa are set to no rank (-1), so we need to track
  // highest rank at each step rather than simply check the final result's
  // rank

  ETaxRank tRank = max(taxRank[taxId1], taxRank[taxId2]),
           hRank = tRank;

  // Test if taxon 1 is a child of taxon 2
  for (int i = parentTaxId[taxId1]; i > 1; i = parentTaxId[i]) {
    if (highestRank != nullptr && hRank < taxRank[i])
      hRank = taxRank[i];

    if (i == taxId2) {
      if (highestRank != nullptr)
        *highestRank = hRank;
      return taxId2;
    }
  }

  // Test if taxon 2 is a child of taxon 1
  hRank = tRank;
  for (int i = parentTaxId[taxId2]; i > 1; i = parentTaxId[i]) {
    if (highestRank != nullptr && hRank < taxRank[i])
      hRank = taxRank[i];

    if (i == taxId1) {
      if (highestRank != nullptr)
        *highestRank = hRank;
      return taxId1;
    }
  }

  // Search for nearest common taxon
  hRank = tRank;
  for (int i = parentTaxId[taxId1]; i > 1; i = parentTaxId[i]) {
    if (highestRank != nullptr && hRank < taxRank[i])
      hRank = taxRank[i];

    for (int j = parentTaxId[taxId2]; j > 1; j = parentTaxId[j]) {
      if (i == j) {
        if (highestRank != nullptr)
          *highestRank = hRank;
        return i;
      }

      if (taxRank[j] == eRank_Superkingdom)
        break; // Don't care about "cellular organisms"/etc. taxa, stop here
    }

    if (taxRank[i] == eRank_Superkingdom)
      break;
  }

  // Default: return root node
  if (highestRank != nullptr)
    *highestRank = taxRank[1]; // eRank_No_Rank
  return 1;
}


/// Performs filtering and other HGT evaluations for a BLAST query set, and
/// records HGT-matched hit pairs with their score.
/// @param results Handle to BLAST search result set [in]
/// @param scope BLAST scope [in]
// ** TODO Add support for null group 2
void CAlienG2::EvaluateResults(CRef<CSearchResultSet> results, CRef<CScope> scope) {
  static int count = 0;

  for (auto &curResult: *results) {
    cout << "\r  Processing query " << ++count << flush;

    CConstRef<CSeq_id> querySeqId = curResult->GetSeqId();
    CBioseq_Handle queryBioseq = scope->GetBioseqHandle(*querySeqId);
    TSeqPos queryLength = queryBioseq.GetBioseqLength();

    // Get this query's hit set
    CConstRef<CSeq_align_set> seqAlignSet = curResult->GetSeqAlign();

    CRef<CSeq_align> bestAlign[2];
    CBioseq_Handle bestBioseq[2];
    double bestEValue[2]   = {0, 0},
           bestBitscore[2] = {0, 0};
    int bestTaxId[2]       = {0, 0};

    // Iterate through this query's hits
    // They are pre-sorted by e-value in ascending order
    for (auto &curAlign: seqAlignSet->Get()) {
      double bitscore, eValue;

      // Test against bitscore, e-value, and coverage thresholds, where applicable
      curAlign->GetNamedScore(CSeq_align::eScore_BitScore, bitscore);
      if (bitscoreThreshold && bitscore < bitscoreThreshold)
        continue;

      curAlign->GetNamedScore(CSeq_align::eScore_EValue, eValue);
      if (eValueThreshold && eValue > eValueThreshold)
        continue;

      const CSeq_id &subjectSeqId = curAlign->GetSeq_id(1);
      CBioseq_Handle subjectBioseq = scope->GetBioseqHandle(subjectSeqId);

      if (coverageThreshold && ((double)curAlign->GetAlignLength() /
          (double)queryLength) < coverageThreshold) {
        continue;
      }

      // Test if subject is filtered by GI
      if (!filteredGIs.empty() && subjectSeqId.IsGi() &&
          filteredGIs.count(subjectSeqId.GetGi()) != 0)
        continue;

      // Test if subject is filtered by accession number
      if (!filteredAccs.empty()) {
        CSeq_id_Handle idForAcc = sequence::GetId(subjectSeqId, *scope,
            sequence::eGetId_ForceAcc);
        if (idForAcc &&
            (filteredAccs.count(idForAcc.GetSeqId()->GetSeqIdString(false)) != 0 ||
             filteredAccs.count(idForAcc.GetSeqId()->GetSeqIdString(true))  != 0))
          continue;
      }

      // Test if subject is filtered by taxon
      int taxId = 0;
      if ((taxId = scope->GetTaxId(subjectSeqId)) <= 1 ||
          IsTaxonInGroup(taxId, filteredTaxa))
        continue;

      // Test if subject's taxon is in HGT group 1 or 2
      bool inGroup1 = IsTaxonInGroup(taxId, group1Taxa),
           inGroup2 = (!inGroup1 && IsTaxonInGroup(taxId, group2Taxa));
      if (!((inGroup1 && bestBitscore[0] == 0) ||
            (inGroup2 && bestBitscore[1] == 0)))
        continue;

      // Test if search order is being followed, if applicable
      if ((group1First && inGroup2 && bestBitscore[0] == 0) ||
          (group2First && inGroup1 && bestBitscore[1] == 0))
        continue;

      // Test if nearest common ancestor passes tax rank threshold. The worse-
      // scoring group is tested against the better-scoring group, not vice-versa
      if (minAncestorRank > eRank_No_Rank &&
          ((inGroup1 && bestBitscore[1]) || (inGroup2 && bestBitscore[0]))) {
        ETaxRank pRank;
        int pTaxId = inGroup1 ?
                     GetCommonAncestor(taxId, bestTaxId[1], &pRank) :
                     GetCommonAncestor(bestTaxId[0], taxId, &pRank);

        if (pTaxId > 1 && pRank < minAncestorRank)
          continue;
      }

      // Store info for the best-scoring hit for the appropriate HGT group
      bestAlign[(int)(inGroup2)]    = curAlign;
      bestBioseq[(int)(inGroup2)]   = subjectBioseq;
      bestEValue[(int)(inGroup2)]   = eValue;
      bestBitscore[(int)(inGroup2)] = bitscore;
      bestTaxId[(int)(inGroup2)]    = taxId;

      // Exit loop if we found all the hits we were looking for
      if (bestBitscore[0] && bestBitscore[1])
        break;
    }
    if (!(bestBitscore[0] && bestBitscore[1]))
      continue; // No hits in one or both HGT groups for this query

    // Test against [group 1 score]:[group 2 score] ratio thresholds, if applicable
    double scoreRatio = bestBitscore[0] / bestBitscore[1];
    if ((scoreRatioMin && scoreRatio < scoreRatioMin) ||
        (scoreRatioMax && scoreRatio > scoreRatioMax))
      continue;

    // HGT index:   [best bitscore for 2nd group] - [best bitscore for 1st group]
    // Alien index: log10([best e-value for 2nd group + 1e-200] /
    //                    [best e-value for 1st group - 1e-200])
    double alienIndex = bestEValue[0] == 0 ?
                        0 : log10((bestEValue[1]) / (bestEValue[0])),
           hgtIndex   = bestBitscore[1] - bestBitscore[0];

    // Test against alien index and HGT index thresholds, if applicable
    if ((alienIndexMin && alienIndex < alienIndexMin) ||
        (alienIndexMax && alienIndex > alienIndexMax))
      continue;

    if ((hgtIndexMin && hgtIndex < hgtIndexMin) ||
        (hgtIndexMax && hgtIndex > hgtIndexMax))
      continue;

    // Store the resultant HGT score, etc. for this query
    hgtData.emplace_back(bestAlign[0], bestAlign[1],
                         sequence::GetId(bestAlign[0]->GetSeq_id(1),
                                         *scope, sequence::eGetId_ForceAcc),
                         sequence::GetId(bestAlign[1]->GetSeq_id(1),
                                         *scope, sequence::eGetId_ForceAcc),
                         bestBioseq[0], bestBioseq[1],
                         alienIndex, hgtIndex, scoreRatio,
                         bestTaxId[0], bestTaxId[1]);
  }
}


/// Write out HGT report file.
/// @param out Output stream to write to
void CAlienG2::CreateReport(CNcbiOstream &out) {
  // Print header
  out << "query_full_title"      << "\t"
      << "first_hit_full_title"  << "\t" << "second_hit_full_title"  << "\t"
      << "first_hit_accession"   << "\t" << "second_hit_accession"   << "\t"
      << "hit1/hit2_score_ratio" << "\t"
      << "hgt_index"             << "\t"
      << "alien_index"           << "\t"
      << "first_hit_len"         << "\t" << "second_hit_len"         << "\t"
      << "first_hit_align_len"   << "\t" << "second_hit_align_len"   << "\t"
      << "first_hit_identities"  << "\t" << "second_hit_identities"  << "\t"
      << "first_hit_coverage"    << "\t" << "second_hit_coverage"    << "\t"
      << "first_hit_range_len"   << "\t" << "second_hit_range_len"   << "\t"
      << "first_hit_range_ratio" << "\t" << "second_hit_range_ratio" << "\t"
      << "first_hit_evalue"      << "\t" << "second_hit_evalue"      << "\t"
      << "first_hit_bitscore"    << "\t" << "second_hit_bitscore"    << "\t"
      << "first_hit_tax_id"      << "\t" << "second_hit_tax_id"      << "\t"
      << "first_hit_species"     << "\t" << "second_hit_species"     << "\t"
      << "first_hit_lineage"     << "\t" << "second_hit_lineage"     << endl;

  for (auto &curHgtDataSet: hgtData) {
    // Unpack entry's data
    CRef<CSeq_align> hit[2];
    CSeq_id_Handle id[2];
    CBioseq_Handle bioseq[2], queryBioseq;
    double alienIndex, hgtIndex, scoreRatio;
    int taxId[2];
    tie(hit[0], hit[1],
        id[0], id[1],
        bioseq[0], bioseq[1],
        alienIndex, hgtIndex, scoreRatio,
        taxId[0], taxId[1]) = curHgtDataSet;

    const CSeq_id &queryId = hit[0]->GetSeq_id(0);
    queryBioseq = bioseq[0].GetScope().GetBioseqHandle(queryId);
    string queryTitle = sequence::GetTitle(queryBioseq); // ** FIXME Deprecated

    // Fetch accessions, titles, scores, and taxonomic lineages
    string idStr[2], title[2];
    TSeqPos length[2];
    double coverage[2], identities[2], evalue[2], bitscore[2];
    TSeqRange seqRange[2];
    vector<int> lineage[2];
    for (int i = 0; i < 2; ++i) {
      if (id[i])
        idStr[i] = id[i].GetSeqId()->GetSeqIdString(true);
      if (idStr[i].empty())
        idStr[i] = hit[i]->GetSeq_id(1).AsFastaString();

      length[i] = bioseq[i].GetBioseqLength();
      coverage[i] = (double)hit[i]->GetAlignLength() /
                    (double)queryBioseq.GetBioseqLength();
      title[i] = sequence::GetTitle(bioseq[i]); // ** FIXME Deprecated
      hit[i]->GetNamedScore(CSeq_align::eScore_IdentityCount, identities[i]);
      hit[i]->GetNamedScore(CSeq_align::eScore_EValue, evalue[i]);
      hit[i]->GetNamedScore(CSeq_align::eScore_BitScore, bitscore[i]);
      seqRange[i] = hit[i]->GetSeqRange(1);

      for (int j = taxId[i]; j > 1; j = parentTaxId[j]) {
        lineage[i].push_back(j);
        if (taxRank[j] == eRank_Superkingdom)
          break; // Don't care about "cellular organisms"/etc. taxa, stop here
      }
    }

    // Print this entry
    out << queryTitle                 << "\t"
        << title[0]                   << "\t" << title[1]                   << "\t"
        << idStr[0]                   << "\t" << idStr[1]                   << "\t"
        << scoreRatio                 << "\t"
        << hgtIndex                   << "\t"
        << alienIndex                 << "\t"
        << length[0]                  << "\t" << length[1]                  << "\t"
        << hit[0]->GetAlignLength()   << "\t" << hit[1]->GetAlignLength()   << "\t"
        << identities[0]              << "\t" << identities[1]              << "\t"
        << coverage[0]                << "\t" << coverage[1]                << "\t"
        << seqRange[0].GetLength()    << "\t" << seqRange[1].GetLength()    << "\t"
        << hit[0]->AlignLengthRatio() << "\t" << hit[1]->AlignLengthRatio() << "\t"
        << evalue[0]                  << "\t" << evalue[1]                  << "\t"
        << bitscore[0]                << "\t" << bitscore[1]                << "\t"
        << taxId[0]                   << "\t" << taxId[1]                   << "\t"
        << sciName[taxId[0]]          << "\t" << sciName[taxId[1]]          << "\t";
    for (int i = 0; i < 2; ++i) {
      out << "[";
      for (auto j = lineage[i].rbegin(); j != lineage[i].rend(); ++j)
        out << "\'"
            << (sciName[*j].empty() ? to_string(*j) : sciName[*j]) << "\'"
            << (*j != taxId[i] ? ", " : "]\t");
    }
    out << endl;        
  }
}


/////////////////////////////////////////////////////////////////////////////
//  Run


int CAlienG2::Run() {
  SetDiagPostLevel(eDiag_Warning);
  SetDiagPostPrefix("AlienG2");

  // Get arguments/parameters
  const CArgs &args = GetArgs();

  bool loadReport = args["task"].AsString().empty() ||
                    args["task"].AsString() == "report";
  EProgram program = !loadReport ? ProgramNameToEnum(args["task"].AsString()) :
                                   (EProgram)0;

#ifndef NCBI_NO_THREADS
  int numBlastThreads = (program == eRPSBlast || program == eRPSTblastn) ? 0 :
                        args["num_threads"].AsInteger();
#else
  // No threads can be set in NON-MT mode
  numBlastThreads = CThreadable::kMinNumThreads;
#endif

  bitscoreThreshold = args["bitscore"].AsDouble();
  eValueThreshold   = args["evalue"].AsDouble();
  coverageThreshold = args["coverage"].AsDouble();
  alienIndexMin     = args["alien_min"].AsDouble();
  alienIndexMax     = args["alien_max"].AsDouble();
  hgtIndexMin       = args["hgt_min"].AsDouble();
  hgtIndexMax       = args["hgt_max"].AsDouble();
  scoreRatioMin     = args["ratio_min"].AsDouble();
  scoreRatioMax     = args["ratio_max"].AsDouble();
  minAncestorRank   = RankNameToEnum(args["ancestor_rank"].AsString());
  group1First       = args["search_order"].AsString() == "group1";
  group2First       = args["search_order"].AsString() == "group2";

  chrono::high_resolution_clock::time_point
      lastTime = chrono::high_resolution_clock::now(),
      curTime;
  chrono::duration<double> timeSpan;

  // Process taxonomy information from taxdump
  cout << "Loading taxdump... " << flush;
  if (!ProcessTaxonomyDump())
    throw std::runtime_error("Could not open taxdump");

  curTime = chrono::high_resolution_clock::now();
  timeSpan = chrono::duration_cast<chrono::duration<double>>(curTime - lastTime);
  lastTime = curTime;
  cout << "done in " << timeSpan.count() << " sec" << endl;

  // Initialize search filters and groups
  SetupFiltersAndGroups();
  if (group1Taxa.empty())
    throw std::runtime_error("HGT group 1 cannot be empty");

  CRef<CObjectManager> objMgr;
  CRef<CScope> lclScope;

  CNcbiOstream &outFile = args["out"].AsOutputFile();

  if (loadReport) {
    cout << "Initializing BLAST report... " << flush;
    // Read input file as BLAST report
    rmtBlast.Reset(new CRemoteBlast(args["in"].AsInputFile()));
    cout << "done" << endl;

    cout << "Analyzing results for HGT candidates..." << endl;
    // Evaluate results
    while (rmtBlast->LoadFromArchive()) { // Load next chunk from file
      CRef<CBlastOptionsHandle> opts = rmtBlast->GetSearchOptions();
      if (opts->GetEvalueThreshold() <= eValueThreshold)
        eValueThreshold = 0; // No need to check e-value threshold

      if (!rmtBlast->IsDbSearch()) // Bl2seq mode
        throw std::runtime_error("Database-free search is not supported");

      // Initialize queries
      CRef<CBlastQueryVector> queries =
        x_ExtractQueries(Blast_QueryIsProtein(opts->GetOptions().GetProgramType()));
      CRef<CScope> scope = queries->GetScope(0);
      _ASSERT(queries);

      // Initialize database
      CRef<CBlast4_database> db = rmtBlast->GetDatabases();
      _ASSERT(db);
      CRef<CSearchDatabase> searchDB(new CSearchDatabase(db->GetName(),
        db->IsProtein() ? CSearchDatabase::eBlastDbIsProtein :
                          CSearchDatabase::eBlastDbIsNucleotide));
      CRef<CBlastDatabaseArgs> dbArgs(new CBlastDatabaseArgs());
      dbArgs->SetSearchDatabase(searchDB);
      CRef<CLocalDbAdapter> dbAdapter;
      InitializeSubject(dbArgs, opts, false, dbAdapter, scope);

      // Evaluate each query's results for HGT candidates
      EvaluateResults(rmtBlast->GetResultSet(), scope);
    }
    curTime = chrono::high_resolution_clock::now();
    timeSpan = chrono::duration_cast<chrono::duration<double>>(curTime - lastTime);
    lastTime = curTime;
    cout << endl << "  Done in " << timeSpan.count() << " sec" << endl;
  }
  else {
    // ** TODO Implement query splitting/batch size
    // Initialize and validate BLAST options
    CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));
    ProcessBlastOptions(opts);
    eValueThreshold = 0;
    opts->Validate();  // Can throw CBlastException::eInvalidOptions
    //opts->GetOptions().DebugDumpText(cout, "opts", 1);
    bool dbIsAA = Blast_SubjectIsProtein(opts->GetOptions().GetProgramType());

    objMgr.Reset(CObjectManager::GetInstance());
    if (!objMgr)
      throw std::runtime_error("Could not initialize object manager");
    lclScope.Reset(new CScope(*objMgr));

    CRef<CSearchDatabase> targetDB(new CSearchDatabase(args["db"].AsString(),
      dbIsAA ? CSearchDatabase::eBlastDbIsProtein :
               CSearchDatabase::eBlastDbIsNucleotide));
    CRef<CBlastDatabaseArgs> dbArgs(new CBlastDatabaseArgs());
    dbArgs->SetSearchDatabase(targetDB);
    CRef<CLocalDbAdapter> dbAdapter;
    InitializeSubject(dbArgs, opts, false, dbAdapter, lclScope);

    cout << "Running " << args["task"].AsString() << "... " << flush;

    // Read input file as FASTA-formatted query file
    SDataLoaderConfig dlconfig(
      Blast_QueryIsProtein(opts->GetOptions().GetProgramType()));
    CBlastInputSourceConfig iconfig(dlconfig);
    CBlastFastaInputSource fastaInput(args["in"].AsInputFile(), iconfig);
    CBlastInput blastInput(&fastaInput);

    // Initialize queries
    TSeqLocVector queryLoc = blastInput.GetAllSeqLocs(*lclScope);
    CRef<IQueryFactory> queryFactory(new CObjMgr_QueryFactory(queryLoc));

    // Execute local BLAST search
    lclBlast.Reset(new CLocalBlast(queryFactory, opts, *targetDB));
    if (numBlastThreads)
      lclBlast->SetNumberOfThreads(numBlastThreads);
    CRef<CSearchResultSet> results = lclBlast->Run();

    curTime = chrono::high_resolution_clock::now();
    timeSpan = chrono::duration_cast<chrono::duration<double>>(curTime - lastTime);
    lastTime = curTime;
    cout << "done in " << timeSpan.count() << " sec" << endl;

    // Get warning messages
    for (auto &curResult: *results) {
      TQueryMessages messages = curResult->GetErrors(eBlastSevWarning);
      if (messages.size() > 0) {
        CConstRef<CSeq_id> querySeqId = curResult->GetSeqId();

        cout << "Query ID: " << (querySeqId.NotEmpty() ? querySeqId->AsFastaString()
                                                       : "Unknown") << endl;
        for (auto &message: messages)
          cout << message->GetMessage() << endl;
      }
    }

    cout << "Analyzing results for HGT candidates..." << endl;
    // Evaluate each query's results for HGT candidates
    EvaluateResults(results, lclScope);

    curTime = chrono::high_resolution_clock::now();
    timeSpan = chrono::duration_cast<chrono::duration<double>>(curTime - lastTime);
    lastTime = curTime;
    cout << endl << "  Done in " << timeSpan.count() << " sec" << endl;
  }

  // Output HGT report
  if (!hgtData.empty()) {
    CreateReport(outFile);
    cout << "HGT report created" << endl;
  }
  else
    cout << "No HGT candidates found" << endl;

  return 0;
}


/////////////////////////////////////////////////////////////////////////////
//  Cleanup


void CAlienG2::Exit() {
  // ...
}


/////////////////////////////////////////////////////////////////////////////
//  MAIN


#ifndef SKIP_DOXYGEN_PROCESSING
int main(int argc, const char* argv[]) {
  // Execute main application function
  return CAlienG2().AppMain(argc, argv);
}
#endif /* SKIP_DOXYGEN_PROCESSING */
