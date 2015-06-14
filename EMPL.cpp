// CS124/224 Final Project - Haplotype Inference using
// Expectation-Maximization Partition-Ligation algorithm.
// By Zeeshan Pirzada


#include <string>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;


// Useful structs 

struct haplotype_pair
{
	haplotype_pair(int a, int b)
	{
		h1 = a;
		h2 = b;
	}

	int h1;
	int h2;
};


// Helper functions

bool no_duplicates(vector<haplotype_pair> pairs, haplotype_pair p)
{
	for (int i = 0; i < pairs.size(); i++)
	{
		haplotype_pair a = pairs[i];

		if ((a.h1 == p.h1 && a.h2 == p.h2)
			|| (a.h1 == p.h2 && a.h2 == p.h1))
			return false;
	}
	return true;
}


int delt(int a, int b)
{
	return (a == b) ? 1 : 2;
}


string haplotypeIntToString(int a, int width)
{
	string ret;

	for (int i = 0; i < width; i++)
	{
		if (a % 2 == 0)
			ret.append("0");
		else
			ret.append("1");

		a = a/2;
	}

	reverse(ret.begin(), ret.end());
	return ret;


}

int haplotypeStringToInt(string hap)
{
	int base = 0;
	int sum = 0;
	
	for (int i = hap.length() - 1; i >= 0 ; i--)
	{
		char c = hap[i];
		int digit = c - '0';
		sum += digit * pow(2, base);
		base++;
	}

	return sum;
}

bool converged(vector<double> old, vector<double> next, double eps)
{
	double diff = 0;
	for (int i = 0; i < old.size(); i++)
	{
		if (old[i] > next[i])
			diff += old[i] - next[i];
		else
			diff += next[i] - old[i];
		cout << diff << endl;
	}


	return (diff <= eps);
}

bool explains(string h1, string h2, string g)
{
	for (int i = 0; i < g.length(); i++)
	{
		if (g[i] == 0)
		{
			if (h1[i] != 0 || h2[i] != 0)
				return false;
		}

		if (g[i] == 1)
		{
			if (h1[i] != 1 || h2[i] != 1)
				return false;
		}

		if (g[i] == 2) // ambiguous
		{
			// g_i = 2 so the haplotype has to be 01 or 10

			if (!(h1[i] ^ h2[i]))
				return false;
		}
	}

	return true;

}




/*
	Return possible haplotype pairs that could
	explain genotype g. This is where the Biology
	is happening.

	To generate the pairs, parse through the genotype.
	At any 0 or 1, this is unambigious so add it to every
	member of the haplotype_pairs. If we get a 2, then for
	every pair in the list, there are 2 choices, 10 or 01. 

	We don't want to have duplicates (or backwards pairs), 
	so we always check for duplicates.

	Note: haplotypes are stored as ints so that indexes can
	be used to store them.
*/


vector<haplotype_pair> getPossibleHapPairs(string g)
{
	vector<haplotype_pair> ret;

	if (g.length() == 0)
		return ret;


	// start with first elem of g.

	if (g[0] == '0' || g[0] == '1')
	{
		
		haplotype_pair hp((g[0] - '0'), (g[0] - '0'));
		ret.push_back(hp);
	}
	else
	{
		haplotype_pair hp1(0, 1);
		haplotype_pair hp2(1, 0);

		ret.push_back(hp1);
		ret.push_back(hp2);

	}


	for (int i = 1; i < g.length(); i++)
	{
		if (g[i] == '0' || g[i] == '1')
		{

			// update existing pairs
			

			for (int j = 0; j < ret.size(); j++)
			{
				// add g[i] to both members of p.


				haplotype_pair p = ret[j];
				
				p.h1 = 2 * p.h1 + (g[i] - '0');
				p.h2 = 2 * p.h2 + (g[i] - '0');
				

				ret[j] = p;

			}
		}
		else  // branch out new pairs
		{
			vector<haplotype_pair> new_pairs;



			for (int j = 0; j < ret.size(); j++)
			{
				// for each pair in ret, make two new
				// pairs representing both branches

				haplotype_pair p1 = ret[j];
				haplotype_pair p2 = ret[j];

				p1.h1 = 2 * p1.h1 + 0;
				p1.h2 = 2 * p1.h2 + 1;

				p2.h1 = 2 * p2.h1 + 1;
				p2.h2 = 2 * p2.h2 + 0;

				new_pairs.push_back(p1);
				new_pairs.push_back(p2);

			}
			ret.clear();
			for (int k = 0; k < new_pairs.size(); k++)
			{

				if (no_duplicates(ret, new_pairs[k]))
					ret.push_back(new_pairs[k]);
			}
		}
	}

	return ret;
}

/*
 	Uses expectation maximization to determine haplotypes.
 	Initially, set the probability of each haplotype to be the same.
 	This could be changed to a random or greedy method later. In the
 	expectation step, go through each genotype and identify all possible
 	haplotype pairs, and compute the probabilities of each of these pairs
 	from the hap_probs vector. Normalize these probabilities and then
 	add to the expected values of all paired haplotypes.

 	In the maximization step, simply rescale the updated probabilities 
 	of the haplotypes.
*/

vector<string> haplotyper_EM(vector<string> genotypes)
{
	
	int len = genotypes[0].length();

	int num_haplotypes = pow(2, len);

	// initialize haplotype probabilities to 1/num_haplotypes

	vector<double> hap_probs;   // hap_probs[i] represents frequency of haplotype i (in binary)

	double p_0 = 1.0/num_haplotypes;

	for (int i = 0; i < num_haplotypes; i++)
	{
		hap_probs.push_back(p_0);	
	}

 

	// main loop

	int num_iters = 0;
	do
	{
		vector<double> hap_expected_values;

		for (int i = 0; i < num_haplotypes; i++)
		{
			hap_expected_values.push_back(0);
		}

		for (int i = 0; i < genotypes.size(); i++)
		{
			
			// Expectation step performed for each genotype

			vector<haplotype_pair> pairs = getPossibleHapPairs(genotypes[i]);
			
			vector<double> pair_probs;
			double sum = 0.0;
			
			// Compute probability of each pair from current haplotype frequencies 
			
		
			for (int j = 0; j < pairs.size(); j++)
			{
				
				int h1 = pairs[j].h1;
				int h2 = pairs[j].h2;

				double p = hap_probs[h1] * hap_probs[h2] * delt(h1, h2);
				pair_probs.push_back(p);
				sum += p;
			}
			

			// add the normalized probabilities to the expected_values vector

			for (int k = 0; k < pairs.size(); k++)
			{
				int h1 = pairs[k].h1;
				int h2 = pairs[k].h2;

				hap_expected_values[h1] += pair_probs[k]/sum;
				hap_expected_values[h2] += pair_probs[k]/sum;

			}

		}

		// Maximization step, scale the hap_probs

		double s = 0.0;
		for (int x = 0; x < hap_expected_values.size(); x++)
		{
			s += hap_expected_values[x];
		}
		for (int y = 0; y < hap_probs.size(); y++)
		{
			hap_probs[y] = hap_expected_values[y]/s;
		}


		num_iters++;

	} while (num_iters < 100);

	// get haplotypes:

	vector<string> ret;

	for (int i = 0; i < hap_probs.size(); i++)
	{
		if (hap_probs[i] > .00001)
		{
			ret.push_back(haplotypeIntToString(i, len));
		}
	}
	return ret;

} 

/*
	function that ligates two blocks together. We do this
	by forming all combinations of the haplotypes in b1 with b2,
	and then running the EM algorithm on the combined data sets.
*/

vector<string> ligateBlocks(vector<string> b1_haplotypes, vector<string> b2_haplotypes, vector<string> genotypes)
{

	// take the cross product of b1 and b2

	vector<string> combined_haps;

	for (int i = 0; i < b1_haplotypes.size(); i++)
	{
		for (int j = 0; j < b2_haplotypes.size(); j++)
		{
			string c = b1_haplotypes[i] + b2_haplotypes[j];
			combined_haps.push_back(c);
		}
	}


	 /*
		now we have a set of haplotypes that could potentially explain
	 	the genotypes, so we run the EM algorithm on them. The implementation
	 	is going to be slightly different here because we are no longer considering
	 	all possible haplotypes like we did in the atomistic units. The main 
	 	difference is we can't index vectors by haplotype. Rather, just use parallel
	 	vectors. */
	 

	 vector<double> hap_probs;
	

	 // initialize probabilities.

	 double p_0 = 1.0/combined_haps.size();
	
	 for (int i = 0; i < combined_haps.size(); i++) 
	 {
	 	hap_probs.push_back(p_0);
	 }

	 cout << p_0 << endl;
	// main loop

	int num_iters = 0;
	do
	{
		vector<double> hap_expected_values;

		for (int i = 0; i < hap_probs.size(); i++)
		{
			hap_expected_values.push_back(0);
		}

		for (int i = 0; i < genotypes.size(); i++)
		{
			
			// Expectation step performed for each genotype
			
			// Compute probability of each haplotype based on genotypes
			
			for (int j = 0; j < combined_haps.size(); j++)
			{
				for (int k = j; k < combined_haps.size(); k++)
				{
					if (explains(combined_haps[j], combined_haps[k], genotypes[i]))
					{
						double to_add = hap_probs[j] * hap_probs[k] * delt(j, k);
						hap_expected_values[j] += to_add;

		//				cout << hap_expected_values[j] << endl;
		//				cout << "Pairs: " << j << " " << k << endl;
					}
				}
			}

		}

		// Maximization step, scale the hap_probs

		double s = 0.0;
		for (int x = 0; x < hap_expected_values.size(); x++)
		{
			s += hap_expected_values[x];
		}
		for (int y = 0; y < hap_probs.size(); y++)
		{
			hap_probs[y] = hap_expected_values[y]/s;
		}

		num_iters++;

	} while (num_iters < 5);

	// get haplotypes:

	vector<string> ret;

	for (int i = 0; i < hap_probs.size(); i++)
	{
		if (hap_probs[i] > .00001)
		{
			ret.push_back(combined_haps[i]);
		}
	}
	return ret;


}


/*
	Helps split the genotypes up for partition-ligation.
*/	

vector<vector<string> > splitGenotypes(vector<string> genotypes, int regionSize)
{

	// split genotypes into regions of 10 SNPs

	int len = genotypes[0].size();

	vector<vector<string> > split_genotypes;

	for (int i = 0; i < len; i += regionSize)
	{
		vector<string> temp;
		for (int j = 0; j < genotypes.size(); j++)
		{
			if (len - i > regionSize)
				temp.push_back(genotypes[j].substr(i, regionSize));
			else
				temp.push_back(genotypes[j].substr(i));
		}
		
		split_genotypes.push_back(temp);
	}
	return split_genotypes;
}


/*
 	
 	Uses expectimization-maximization and partition ligation.
	The EM algorithm is O(2^k), so it sucks for K>15ish. Partition
	ligation overcomes this by splitting the genotypes into sub problems
	and running EM on the subproblems and then combining the
	solutions.

	When we ligate two blocks, assume with N1 and N2 haplotypes phased
	from each region. In order to combine the two blocks, we have to 
	take cross product to get all possibilities. The frequencies will
	multiply, and so the threshhold should be reduced by a factor of 100?

*/

vector<string> haplotyper_EMPL(vector<string> genotypes)
{

	int len = genotypes[0].length();

	// if small number of SNPs, no point using PL

	if (len < 15)
		return haplotyper_EM(genotypes);
	
	// split the genotypes into blocks of at 10 and solve them.

	vector<vector<string> > split_genotypes = splitGenotypes(genotypes, 5);

	vector<vector<string> > hap_results;

	for (int i = 0; i < split_genotypes.size(); i++)
	{
		hap_results.push_back(haplotyper_EM(split_genotypes[i]));
	}

	// recursively apply ligate function in divide and conquer fashion

	for (;;)
	{
		if (hap_results.size() == 1)
			break;

		vector<vector<string> > partial_hap_results;
		vector<vector<string> > partial_genotypes;

		for (int i = 0; i < hap_results.size(); i += 2)
		{

			// recombine genotype blocks i and i+1.
			int n = split_genotypes[0].size();
			vector<string> combined_genotypes;
			
			for (int j = 0; j < n; j++)
			{
				combined_genotypes.push_back( split_genotypes[i][j] + split_genotypes[i+1][j] );			
				
			}
			partial_genotypes.push_back(combined_genotypes);
			
			// ligate the blocks.
			partial_hap_results.push_back(ligateBlocks(hap_results[i], hap_results[i+1], partial_genotypes[i/2]));

		}
		
		split_genotypes = partial_genotypes;
		hap_results = partial_hap_results;
	}

	
	return hap_results[0];


}

// generate genotypes

vector<string> generateGenotypes(int numGenotypes, int numSNPS)
{
	vector<string> ret;

	for (int i = 0; i < numGenotypes; i++)
	{
		string a;
		for (int j = 0; j < numSNPS; j++)
		{
			int pos = rand() % 3;
			a.append(to_string(pos));
		}
		ret.push_back(a);
	}
	return ret;
}


int main(int argc, char *argv[])
{
	
	srand (time (NULL));

	if (argc != 3)
	{
		cout << "Usage: ./empl <num-people> <num-SNPs>" << endl;
		return 0;
	}

	int num_people = atoi(argv[1]);
	int num_SNPs = atoi(argv[2]);

	vector<string> g = generateGenotypes(num_people, num_SNPs);

	vector<string> h = haplotyper_EMPL(g);
	cout << "Size of solution from EM: " << h.size() << endl;

	vector<string> hh = haplotyper_EM(g);
	cout << "Size of solution from EM-PL: " << hh.size() << endl;


}



