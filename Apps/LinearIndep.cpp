#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Element.h"
#include "LRSpline/Meshline.h"

using namespace Go;
using namespace LR;
using namespace std;

int main(int argc, char **argv) {
#ifdef TIME_LRSPLINE
	Profiler prof(argv[0]);
#endif

	// set default parameter values
	bool verbose                  = false;
	bool dumpFile                 = false;
	bool one_by_one               = false;
	double beta                   = 0.10;
	enum refinementStrategy strat = LR_SAFE;
	int mult                      = 1;
	int symmetry                  = 1;
	char *inputFileName           = NULL;
	char *refineFileName          = NULL;
	string parameters(" parameters: \n" \
	                  "   -refine <s>   test refinement given by input file <s>\n"\
	                  "   -verbose      verbose output\n"\
	                  "   -dumpfile     writes an eps- and lr-file of the LR-mesh\n"\
	                  "   -one-by-one   does the meshline insertions one at a time, checking for independence at every step\n"\
	                  "   -help         display (this) help screen\n"
	                  " Refinement file syntax\n"\
	                  "   <beta> <multiplicity> <strategy> <symmetry>\n"\
	                  "   <number of elements to follow>\n"\
	                  "   <highest error element index #0>\n"\
	                  "   <highest error element index #1>\n"\
	                  "   <highest error element index #2>\n"\
	                  "                  ...              \n");
	
	// read input
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i], "-dumpfile") == 0)
			dumpFile = true;
		else if(strcmp(argv[i], "-verbose") == 0)
			verbose = true;
		else if(strcmp(argv[i], "-one-by-one") == 0)
			one_by_one = true;
		else if(strcmp(argv[i], "-refine") == 0)
			refineFileName = argv[++i];
		else if(strcmp(argv[i], "-help") == 0) {
			cout << "usage: " << argv[0] << "[parameters] <refine inputfile>" << endl << parameters;
			exit(0);
		} else {
			if(inputFileName != NULL) {
				cerr << "usage: " << argv[0] << endl << parameters;
				exit(1);
			} else {
				inputFileName = argv[i];
			}
		}
	}

	// do some error testing on input
	if(inputFileName == NULL) {
		cerr << "ERROR: Specify input file name\n";
		cerr << "usage: " << argv[0] << "[parameters] <refine inputfile>" << endl << parameters;
		exit(3);
	}

	LRSplineSurface *lr;

	ifstream inputfile;
	inputfile.open(inputFileName);
	if(!inputfile.is_open()) {
		cerr << "Error: could not open file " << inputFileName << endl;
		exit(3);
	}
	lr = new LRSplineSurface();
	inputfile >> *lr;

	bool isLinearIndep = true;
	int meshlineFail   = -1;
	if(refineFileName != NULL) {
		ifstream refineFile;
		refineFile.open(refineFileName);
		if(!refineFile.is_open()) {
			cerr << "Error: could not open file " << refineFileName << endl;
			exit(3);
		}
		cout << "Reading refine file " << refineFileName << "...\n";
		int strat_index;
		vector<int> sorted_list;
		refineFile >> beta;
		refineFile >> mult >> strat_index >> symmetry;
		if(strat_index == 0) 
			strat = LR_SAFE;
		else if(strat_index == 1) 
			strat = LR_MINSPAN;
		else if(strat_index == 2) 
			strat = LR_ISOTROPIC_EL;
		else if(strat_index == 3) 
			strat = LR_ISOTROPIC_FUNC;
		else {
			cerr << "Error: Invalid refinement strategy " << strat_index << " (valid input: 0,1,2,3)\n";
			exit(4);
		}
		int a;
		int N;
		ws(refineFile);
		refineFile >> N;
		for(int i=0; i<N; i++) {
			refineFile >> a;
			ws(refineFile);
			sorted_list.push_back(a);
		}
		cout << "done\n";

		cout << "calling LRSplineSurface::refine(...)\n";

		LRSplineSurface *lr_original = NULL;
		if(one_by_one)
			lr_original = lr->copy();
			
		cout << setprecision(16);
		vector<Meshline*> *newLines = new vector<Meshline*>();
		lr->refine(sorted_list, beta, mult, strat, symmetry, newLines);
		cout << "Number of new mesh lines: " << newLines->size() << endl;
		for(unsigned int i=0; i<newLines->size(); i++) {
			newLines->at(i)->writeMore(cout);
			cout << endl;
		}
		
		if(one_by_one) {
			cout << "Inserting meshlines one by one...\n";
			for(unsigned int i=0; i<newLines->size(); i++) {
				Meshline *m = newLines->at(i);
				printf("Nelements = %5d Nbasis = %5d, ", lr_original->nElements(), lr_original->nBasisFunctions());

				if(m->type_ == NEWLINE)
					cout << "Inserting line: " << *m << endl;
				else if(m->type_ == ELONGATION)
					cout << "Expanding line: " << *m << endl;
				else if(m->type_ == MERGING)
					cout << "Merging lines : " << *m << endl;

				if(m->is_spanning_u())
					lr_original->insert_const_v_edge(m->const_par_, m->start_, m->stop_, m->multiplicity_);
				else
					lr_original->insert_const_u_edge(m->const_par_, m->start_, m->stop_, m->multiplicity_);
				if( ! lr_original->isLinearIndepByMappingMatrix(false) ) {
					printf("Nelements = %5d Nbasis = %5d \n", lr_original->nElements(), lr_original->nBasisFunctions());
					lr = lr_original;
					isLinearIndep = false;
					break;
				}
			}
		}
	}

	if(!one_by_one)
		isLinearIndep = lr->isLinearIndepByMappingMatrix(verbose);

	if(dumpFile) {
		ofstream meshfile;
		meshfile.open("mesh.eps");
		lr->writePostscriptMesh(meshfile);
		meshfile.close();

		ofstream functionfile;
		functionfile.open("functions.eps");
		lr->writePostscriptFunctionSpace(functionfile);
		functionfile.close();

		ofstream domainfile;
		domainfile.open("domain.eps");
		lr->writePostscriptElements(domainfile, 10, 10);
		domainfile.close();

		ofstream controlmesh;
		controlmesh.open("controlmesh.eps");
		lr->writePostscriptMeshWithControlPoints(controlmesh, 10, 10);
		controlmesh.close();
	
		ofstream lrfile;
		lrfile.open("lrspline.txt");
		lrfile << *lr << endl;
		lrfile.close();

		cout << endl;
		cout << "Written mesh to mesh.eps, functions.eps, domain.eps, controlmesh.eps and lrspline.txt\n";
	}

	if(isLinearIndep) {
		cout << "Mesh is linear independent\n";
		exit(0);
	} else {
		cout << "Linear dependent mesh!\n";
		exit(1);
	}

}
