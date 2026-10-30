#include <MPO.hpp>


namespace spinchain{


MPO::MPO(const SpinChain& model)
{
    std::vector<BondTerm> terms = model.bondTerms();
    MatrixXcd onsite = model.onsiteTerm();
    Index d = model.site.dim;
    MatrixXcd identity = MatrixXcd::Identity(d,d);

    Index D = terms.size() + 2;
    Index last = D - 1;
    size_t length = model.chainlength;

    for (size_t i=0; i<length; i++){
        bool isFirst = (i == 0);
        bool isLast = (i + 1 == length);

        Tensor4 site(isFirst ? 1 : D, isLast ? 1 : D, d, d);
        site.setZero();

        // writes an operator to entry (a,b) of the matrix in the class documentation,
        // ignoring the rows and columns which the chain boundaries remove
        auto place = [&](Index a, Index b, const MatrixXcd& op){
            if ((isFirst && a != 0) || (isLast && b != last)){
                return;
            }
            for (Index s=0; s<d; s++){
                for (Index t=0; t<d; t++){
                    site(isFirst ? 0 : a, isLast ? 0 : b, s, t) = op(s,t);
                }
            }
        };

        place(0, 0, identity);
        place(last, last, identity);
        place(0, last, onsite);
        for (size_t k=0; k<terms.size(); k++){
            place(0, k+1, terms[k].left);
            place(k+1, last, terms[k].coupling*terms[k].right);
        }
        W.push_back(site);
    }
}

} //namespace spinchain
