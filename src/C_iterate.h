
// void show(dblIteratorList x);
// void show(std::vector<double> x);

dblIteratorList getStarts(dblList x);
dblIteratorList getEnds(dblList x);

void increase(dblIteratorList& ii, int& changed,
              const dblIteratorList& _starts_, 
              const dblIteratorList& _ends_, 
              int pos = -1);

bool finished(dblIteratorList ii, 
              const dblIteratorList _ends_);

long long int ll(dblList x);
