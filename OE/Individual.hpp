#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <utility>
#include <string>
#include <ostream>
#include <istream>

// The EvoLearn program runs evolutionary simulations
// Copyright (C) 2025  Olof Leimar
// See Readme.md in top repo directory for copyright notice

//************************ Struct Individual *****************************

// This struct represents an individual with genotype of type GenType,
// phenotype of type PhenType

template<typename GenType, typename PhenType>
struct Individual {
// public:
    using gen_type = GenType;
    using gam_type = typename gen_type::gam_type;
    using mut_rec_type = typename gen_type::mut_rec_type;
    using rho_vec_type = typename gen_type::rho_vec_type;
    using phen_type = PhenType;
    Individual() {}
    Individual(int nsd, gen_type&& g, phen_type&& ph) :
        genotype{g},
        phenotype{nsd, g} {}
    // Construct individual from one gamete
    Individual(int nsd, gam_type&& gam) :
        genotype(std::forward<gam_type>(gam)),
        phenotype(nsd, genotype) {}
    Individual(int nsd, const gam_type& gam) :
        genotype(gam),
        phenotype(nsd, genotype) {}
    // Construct individual from maternal and paternal gametes
    Individual(int nsd, gam_type&& mat_gam, gam_type&& pat_gam) :
        genotype(std::forward<gam_type>(mat_gam),
                 std::forward<gam_type>(pat_gam)),
        phenotype(nsd, genotype) {}
    void Assign(gam_type&& gam);
    void Assign(gam_type&& mat_gam, gam_type&& pat_gam);
    void Assign(gen_type&& g, phen_type&& ph);
    gam_type GetGamete(mut_rec_type& mr) const
    { return genotype.GetGamete(mr); }
    gam_type GetGamete(mut_rec_type& mr, const rho_vec_type& rho) const
    { return genotype.GetGamete(mr, rho); }
    const gen_type& Genotype() const { return genotype; }
    gen_type& Genotype() { return genotype; }
    const phen_type& Phenotype() const { return phenotype; }
    phen_type& Phenotype() { return phenotype; }
    int SubPopNum() const { return phenotype.spn; }
    void SetSubPopN(int a_spn) { phenotype.spn = a_spn; }
    bool Alive() const { return phenotype.alive; }
    void SetAlive() { phenotype.alive = true; }
    void SetDead() { phenotype.alive = false; }
    bool Female() const { return phenotype.Female(); }
    void SetFemale(bool female) { phenotype.female = female; }
    static std::string ColHeads(int n_loci, int nsd);
    // public data members
    gen_type genotype;
    phen_type phenotype;
};

// Construct individual from one gamete and a spn
template<typename GenType, typename PhenType>
void Individual<GenType, PhenType>::Assign(gam_type&& gam)
{
    genotype.Assign(gam);
    phenotype.Assign(genotype);
}

// Construct individual by assigning maternal and paternal gametes
template<typename GenType, typename PhenType>
void Individual<GenType, PhenType>::Assign(gam_type&& mat_gam,
    gam_type&& pat_gam)
{
    genotype.Assign(mat_gam, pat_gam);
    phenotype.Assign(genotype);
}

// Construct an individual by assigning all its data
template<typename GenType, typename PhenType>
void Individual<GenType, PhenType>::Assign(gen_type&& g, phen_type&& ph)
{
    genotype = g;
    phenotype = ph;
}


template<typename GenType, typename PhenType>
std::string Individual<GenType, PhenType>::ColHeads(int n_loci,
                                                    int nsd)
{
    std::string col_hds = gen_type::ColHeads(n_loci);
    col_hds += "\t";
    col_hds += phen_type::ColHeads(nsd);
    return col_hds;
}


//---------------------------------------------------------------------
// Ouput and input of Individual<T> objects

template<typename GenType, typename PhenType>
std::ostream& operator<<(std::ostream& ostr,
                         const Individual<GenType, PhenType>& ind)
{
    ostr << ind.Genotype() << '\t'
         << ind.Phenotype();
    return ostr;
}

template<typename GenType, typename PhenType>
std::istream& operator>>(std::istream& istr,
                         Individual<GenType, PhenType>& ind)
{
    istr >> ind.Genotype()
         >> ind.Phenotype();
    return istr;
}


#endif // INDIVIDUAL_HPP
