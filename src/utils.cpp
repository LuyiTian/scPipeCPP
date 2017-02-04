#include "utils.h"

using namespace std;

string join_path(const string p1, const string p2)
{
    auto sep = '/';
	return p1.back() == sep ? p1 + p2 : p1 + sep + p2;
}

int hamming_distance(string A, string B)
{
    int dist = 0;
    for (int i = 0; i < A.length(); ++i)
    {
        if (A[i] != B[i])
        {
            dist++;
        }
    }
    return dist;
}

void check_file_exists(string fn)
{
    ifstream f(fn.c_str());
    if (f.good()) {
        f.close();
    } else {
        f.close();
        throw invalid_argument("cannot open file: "+ fn + "\n");
    }   
}

// tally the element in vector
unordered_map<string, int> vector_counter(vector<string> v)
{
    unordered_map<string, int> counter;
    for(auto const& val: v)
    {
        if (counter.find(val) != counter.end())
        {
            counter[val] ++;
        }
        else
        {
            counter[val] = 1;
        }
    }

    return counter;
}

// Constructors for Interval
Interval::Interval(int s, int e): st(s), en(e), snd(0) {}
Interval::Interval(int s, int e, int sd): st(s), en(e), snd(sd) {}

int Interval::overlap(int st1, int en1)
{
  int s, e;
  if (en1 < st)
  {
    return -1;
  }
  else if (st1 > en)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

bool operator < (const Interval &L, const Interval &R)
{
  return L.en < R.st;
}

bool operator > (const Interval &L, const Interval &R)
{
  return L.st > R.en;
}

vector<string> &split(const string &s, char delim, vector<string> &elems)
{
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim)
{
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}