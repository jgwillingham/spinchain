#include <Solver.hpp>
#include <numeric>


namespace spinchain{


Solver::Solver(const SpinChain& model) : model(model)
{
}



Vector3cd Solver::magnetization() const
{
    auto total = [this](const MatrixXcd& op){
        std::vector<std::complex<double>> values = localExpectations(op);
        return std::accumulate(values.begin(), values.end(), std::complex<double>(0.0));
    };

    Vector3cd magnetization;
    magnetization << total(model.site.Sx), total(model.site.Sy), total(model.site.Sz);
    return magnetization;
}

} //namespace spinchain
