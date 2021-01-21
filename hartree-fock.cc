/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

// standard C++ headers
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <chrono>

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Libint Gaussian integrals library
#include <libint2.hpp>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        Matrix;  // import dense, dynamically sized Matrix type from Eigen;
                 // this is a matrix with row-major storage (http://en.wikipedia.org/wiki/Row-major_order)
                 // to meet the layout of the integrals returned by the Libint integral library

struct Atom {
    int atomic_number;
    double x, y, z;
};

std::vector<libint2::Atom> read_geometry(const std::string& filename);
std::vector<libint2::Shell> make_sto3g_basis(const std::vector<libint2::Atom>& atoms);
std::vector<libint2::Shell> make_cc_pvdz_basis(const std::vector<libint2::Atom>& atoms);
//size_t nbasis(const std::vector<libint2::Shell>& phells);
//### commented def of below
//###std::vector<size_t> map_shell_to_basis_function(const libint2::BasisSet& obs);
//#std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells);
Matrix compute_soad(const std::vector<libint2::Atom>& atoms);

//#### try new declaration?
Matrix compute_1body_ints(const libint2::BasisSet& obs,
                          libint2::Operator obtype,
                          const std::vector<libint2::Atom>& atoms = std::vector<libint2::Atom>());
/*Matrix compute_1body_ints(const libint2::BasisSet& obs,
                          libint2::Operator obtype,
                          const std::vector<Atom>& atoms)*/
/*##Matrix compute_1body_ints(const std::vector<libint2::Shell>& shells,
                          libint2::Operator t,
                          const std::vector<Atom>& atoms = std::vector<Atom>());*/
// simple-to-read, but inefficient Fock builder; computes ~16 times as many ints as possible
Matrix compute_2body_fock_simple(const libint2::BasisSet& obs,
                                 const Matrix& D);
/*##Matrix compute_2body_fock_simple(const std::vector<libint2::Shell>& shells,
                                 const Matrix& D);
								 */
// an efficient Fock builder; *integral-driven* hence computes permutationally-unique ints once
Matrix compute_2body_fock(const libint2::BasisSet& obs,
                                 const Matrix& D);
/*##
Matrix compute_2body_fock(const std::vector<libint2::Shell>& shells,
                                 const Matrix& D);
								 */

int main(int argc, char *argv[]) {

  using std::cout;
  using std::cerr;
  using std::endl;

  //using libint2::Shell;
  using libint2::BasisSet;
  using libint2::Engine;
  using libint2::Operator;
  using libint2::Atom; //### tested this out at some point

  //### duplicate? using libint2::BasisSet;

  try {

    /*** =========================== ***/
    /*** initialize molecule         ***/
    /*** =========================== ***/

    // read geometry from a file; by default read from h2o.xyz, else take filename (.xyz) from the command line
    const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";
    std::vector<libint2::Atom> atoms = read_geometry(filename);
	const auto basisname = (argc > 2) ? argv[2] : "STO-3g"; //aug-cc-pVDZ

    // count the number of electrons
    auto nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i)
      nelectron += atoms[i].atomic_number;
    const auto ndocc = nelectron / 2;

    // compute the nuclear repulsion energy
    auto enuc = 0.0;
    for (auto i = 0; i < atoms.size(); i++)
      for (auto j = i + 1; j < atoms.size(); j++) {
        auto xij = atoms[i].x - atoms[j].x;
        auto yij = atoms[i].y - atoms[j].y;
        auto zij = atoms[i].z - atoms[j].z;
        auto r2 = xij*xij + yij*yij + zij*zij;
        auto r = sqrt(r2);
        enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;
      }
    cout << "\tNuclear repulsion energy = " << enuc << endl;

    /*** =========================== ***/
    /*** create basis set            ***/
    /*** =========================== ***/

      //####auto shells = make_sto3g_basis(atoms); // not necessary?
//    auto shells = make_cc_pvdz_basis(atoms);
    BasisSet obs(basisname, atoms); //TODO: finish here, create basisname, switch shells fro basisset
	//### took out Libint2:: prefix


    size_t nao = 0;
    for (auto s=0; s<obs.size(); ++s) //#### replace by obs
      nao += obs[s].size();

    /*** =========================== ***/
    /*** compute 1-e integrals       ***/
    /*** =========================== ***/

    // initializes the Libint integrals library ... now ready to compute
    libint2::initialize();

    // compute overlap integrals
    auto S = compute_1body_ints(obs, Operator::overlap, atoms); //#####
    //# auto S = compute_1body_ints(shells, Operator::overlap);
    cout << "\n\tOverlap Integrals:\n";
    cout << S << endl;

    // compute kinetic-energy integrals
    auto T = compute_1body_ints(obs, Operator::kinetic, atoms); //#####
    //##auto T = compute_1body_ints(shells, Operator::kinetic);
    cout << "\n\tKinetic-Energy Integrals:\n";
    cout << T << endl;

    // compute nuclear-attraction integrals
    //##Matrix V = compute_1body_ints(shells, Operator::nuclear, atoms);
    Matrix V = compute_1body_ints(obs, Operator::nuclear, atoms); //#####
    cout << "\n\tNuclear Attraction Integrals:\n";
    cout << V << endl;

    // Core Hamiltonian = T + V
    Matrix H = T + V;
    cout << "\n\tCore Hamiltonian:\n";
    cout << H << endl;

    // T and V no longer needed, free up the memory
    T.resize(0,0);
    V.resize(0,0);

    /*** =========================== ***/
    /*** build initial-guess density ***/
    /*** =========================== ***/

    const auto use_hcore_guess = false;  // use core Hamiltonian eigenstates to guess density?
                                         // set to true to match the result of versions 0, 1, and 2 of the code
                                         // HOWEVER !!! even for medium-size molecules hcore will usually fail !!!
                                         // thus set to false to use Superposition-Of-Atomic-Densities (SOAD) guess
    Matrix D;
    if (use_hcore_guess) { // hcore guess
      // solve H C = e S C
      Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H, S);
      auto eps = gen_eig_solver.eigenvalues();
      auto C = gen_eig_solver.eigenvectors();
      cout << "\n\tInitial C Matrix:\n";
      cout << C << endl;

      // compute density, D = C(occ) . C(occ)T
      auto C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();
    }
    else {  // SOAD as the guess density, assumes STO-nG basis
      D = compute_soad(atoms);
    }

    cout << "\n\tInitial Density Matrix:\n";
    cout << D << endl;

    /*** =========================== ***/
    /*** main iterative loop         ***/
    /*** =========================== ***/

    const auto maxiter = 300;
    const auto conv = 1e-12;
    auto iter = 0;
    auto rmsd = 0.0;
    auto ediff = 0.0;
    auto ehf = 0.0;
    do {
      const auto tstart = std::chrono::high_resolution_clock::now();
      ++iter;

      // Save a copy of the energy and the density
      auto ehf_last = ehf;
      auto D_last = D;

      // build a new Fock matrix
      auto F = H;
      //F += compute_2body_fock_simple(shells, D);
      //##F += compute_2body_fock(shells, D);
      F += compute_2body_fock(obs, D);

      if (iter == 1) {
        cout << "\n\tFock Matrix:\n";
        cout << F << endl;
      }

      // solve F C = e S C
      Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
      auto eps = gen_eig_solver.eigenvalues();
      auto C = gen_eig_solver.eigenvectors();

      // compute density, D = C(occ) . C(occ)T
      auto C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();

      // compute HF energy
      ehf = 0.0;
      for (auto i = 0; i < nao; i++)
        for (auto j = 0; j < nao; j++)
          ehf += D(i,j) * (H(i,j) + F(i,j));

      // compute difference with last iteration
      ediff = ehf - ehf_last;
      rmsd = (D - D_last).norm();

      const auto tstop = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> time_elapsed = tstop - tstart;

      if (iter == 1)
        std::cout <<
        "\n\n Iter        E(elec)              E(tot)               Delta(E)             RMS(D)         Time(s)\n";
      printf(" %02d %20.12f %20.12f %20.12f %20.12f %10.5lf\n", iter, ehf, ehf + enuc,
             ediff, rmsd, time_elapsed.count());

    } while (((fabs(ediff) > conv) || (fabs(rmsd) > conv)) && (iter < maxiter));

    printf("** Hartree-Fock energy = %20.12f\n", ehf + enuc);

    libint2::finalize(); // done with libint

  } // end of try block; if any exceptions occurred, report them and exit cleanly

  catch (const char* ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  }
  catch (std::string& ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  }
  catch (std::exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }
  catch (...) {
    cerr << "caught unknown exception\n";
    return 1;
  }

  return 0;
}

// this reads the geometry in the standard xyz format supported by most chemistry software
std::vector<libint2::Atom> read_dotxyz(std::istream& is) {
  // line 1 = # of atoms
  size_t natom;
  is >> natom;
  // read off the rest of line 1 and discard
  std::string rest_of_line;
  std::getline(is, rest_of_line);

  // line 2 = comment (possibly empty)
  std::string comment;
  std::getline(is, comment);

  std::vector<libint2::Atom> atoms(natom);
  for (auto i = 0; i < natom; i++) {
    std::string element_label;
    double x, y, z;
    is >> element_label >> x >> y >> z;

    // .xyz files report element labels, hence convert to atomic numbers
    int Z;
    if (element_label == "H")
      Z = 1;
    else if (element_label == "C")
      Z = 6;
    else if (element_label == "N")
      Z = 7;
    else if (element_label == "O")
      Z = 8;
    else if (element_label == "F")
      Z = 9;
    else if (element_label == "S")
      Z = 16;
    else if (element_label == "Cl")
      Z = 17;
    else {
      std::cerr << "read_dotxyz: element label \"" << element_label << "\" is not recognized" << std::endl;
      throw "Did not recognize element label in .xyz file";
    }

    atoms[i].atomic_number = Z;

    // .xyz files report Cartesian coordinates in angstroms; convert to bohr
    const auto angstrom_to_bohr = 1 / 0.52917721092; // 2010 CODATA value
    atoms[i].x = x * angstrom_to_bohr;
    atoms[i].y = y * angstrom_to_bohr;
    atoms[i].z = z * angstrom_to_bohr;
  }

  return atoms;
}

std::vector<libint2::Atom> read_geometry(const std::string& filename) {

  std::cout << "Will read geometry from " << filename << std::endl;
  std::ifstream is(filename);
  assert(is.good());

  // to prepare for MPI parallelization, we will read the entire file into a string that can be
  // broadcast to everyone, then converted to an std::istringstream object that can be used just like std::ifstream
  std::ostringstream oss;
  oss << is.rdbuf();
  // use ss.str() to get the entire contents of the file as an std::string
  // broadcast
  // then make an std::istringstream in each process
  std::istringstream iss(oss.str());

  // check the extension: if .xyz, assume the standard XYZ format, otherwise throw an exception
  if ( filename.rfind(".xyz") != std::string::npos)
    return read_dotxyz(iss);
  else
    throw "only .xyz files are accepted";
}

/*###
std::vector<libint2::Shell> make_sto3g_basis(const std::vector<Atom>& atoms) {

  using libint2::Shell;

  std::vector<Shell> shells;

  for(auto a=0; a<atoms.size(); ++a) {

    // STO-3G basis set
    // cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical Physics 51, 2657 (1969)
    //       doi: 10.1063/1.1672392
    // obtained from https://bse.pnl.gov/bse/portal
    switch (atoms[a].atomic_number) {
      case 1: // Z=1: hydrogen
        shells.push_back(
            {
              {3.425250910, 0.623913730, 0.168855400}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
        );
        break;

      case 6: // Z=6: carbon
        shells.push_back(
            {
              {71.616837000, 13.045096000, 3.530512200},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {2.941249400, 0.683483100, 0.222289900},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {2.941249400, 0.683483100, 0.222289900},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;

      case 7: // Z=7: nitrogen
        shells.push_back(
            {
              {99.106169000, 18.052312000, 4.885660200},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {3.780455900, 0.878496600, 0.285714400},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
          {3.780455900, 0.878496600, 0.285714400},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;

      case 8: // Z=8: oxygen
        shells.push_back(
            {
              {130.709320000, 23.808861000, 6.443608300},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {5.033151300, 1.169596100, 0.380389000},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {5.033151300, 1.169596100, 0.380389000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;

    case 9: // Z=9: fluorine
      shells.push_back(
          {
            {166.6791300, 30.3608120, 8.2168207},
            {
              {0, false, {0.15432897, 0.53532814, 0.44463454}}
            },
            {{atoms[a].x, atoms[a].y, atoms[a].z}}
          }
      );
      shells.push_back(
          {
            {6.4648032, 1.5022812, 0.4885885},
            {
              {0, false, {-0.09996723, 0.39951283, 0.70011547}}
            },
            {{atoms[a].x, atoms[a].y, atoms[a].z}}
          }
      );
      shells.push_back(
          {
            {6.4648032, 1.5022812, 0.4885885},
            { // contraction 0: p shell (l=1), spherical=false
              {1, false, {0.15591627, 0.60768372, 0.39195739}}
            },
            {{atoms[a].x, atoms[a].y, atoms[a].z}}
          }
      );
      break;

      default:
        throw "do not know STO-3G basis for this Z";
    }

  }

  return shells;
} */

/*##
std::vector<libint2::Shell> make_cc_pvdz_basis(const std::vector<Atom>& atoms) {
    using libint2::Shell;

    std::vector<Shell> shells;

    for(auto a=0; a<atoms.size(); ++a) {

      // cc-pVDZ basis set
      // obtained from https://bse.pnl.gov/bse/portal
      switch (atoms[a].atomic_number) {
        case 1: // Z=1: hydrogen
          shells.push_back(
              {
                {13.0100000, 1.9620000, 0.4446000, 0.1220000}, // exponents of primitive Gaussians
                {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                  {0, false, {0.0196850, 0.1379770, 0.4781480, 0.5012400}}
                },
                {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
              }
          );
          shells.push_back(
              {
                {0.1220000},
                {
                  {0, false, {1.0000000}}
                },
                {{atoms[a].x, atoms[a].y, atoms[a].z}}
              }
          );
          shells.push_back(
              {
                {0.7270000},
                {
                  {1, false, {1.0000000}}
                },
                {{atoms[a].x, atoms[a].y, atoms[a].z}}
              }
          );

          break;

      case 6: // Z=6: carbon
        shells.push_back(
            {
              {6665.0000000, 1000.0000000, 228.0000000, 64.7100000, 21.0600000, 7.4950000, 2.7970000, 0.5215000, 0.1596000}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {0.0006920, 0.0053290, 0.0270770, 0.1017180, 0.2747400, 0.4485640, 0.2850740, 0.0152040, -0.0031910}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
        );
        shells.push_back(
            {
              {6665.0000000, 1000.0000000, 228.0000000, 64.7100000, 21.0600000, 7.4950000, 2.7970000, 0.5215000, 0.1596000}, // exponents of primitive Gaussians
              {
                {0, false, {-0.0001460, -0.0011540, -0.0057250, -0.0233120, -0.0639550, -0.1499810, -0.1272620, 0.5445290, 0.5804960}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.1596000},
              {
                {0, false, {1.0000000}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {9.4390000, 2.0020000, 0.5456000, 0.1517000},
              {
                {1, false, {0.0381090, 0.2094800, 0.5085570, 0.4688420}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.1517000},
              {
                {1, false, {1.0000000}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.5500000},
              {
                {2, true, {1.0000000}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );

          break;

      case 8: // Z=8: oxygen
        shells.push_back(
            {
              {11720.0000000, 1759.0000000, 400.8000000, 113.7000000, 37.0300000, 13.2700000, 5.0250000, 1.0130000, 0.3023000}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {0.0007100, 0.0054700, 0.0278370, 0.1048000, 0.2830620, 0.4487190, 0.2709520, 0.0154580, -0.0025850}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
        );
        shells.push_back(
            {
              {11720.0000000, 1759.0000000, 400.8000000, 113.7000000, 37.0300000, 13.2700000, 5.0250000, 1.0130000, 0.3023000}, // exponents of primitive Gaussians
              {
                {0, false, {-0.0001600, -0.0012630, -0.0062670, -0.0257160, -0.0709240, -0.1654110, -0.1169550, 0.5573680, 0.5727590}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.3023000},
              {
                {0, false, {1.0000000}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {17.7000000, 3.8540000, 1.0460000, 0.2753000},
              {
                {1, false, {0.0430180, 0.2289130, 0.5087280, 0.4605310}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.2753000},
              {
                {1, false, {1.0000000}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.1850000},
              {
                {2, true, {1.0000000}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );

          break;

      case 9: // Z=9: fluorine
        shells.push_back(
            {
              {14710.0000000, 2207.0000000, 502.8000000, 142.6000000, 46.4700000, 16.7000000, 6.3560000, 1.3160000, 0.3897000}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {0.0007210, 0.0055530, 0.0282670, 0.1064440, 0.2868140, 0.4486410, 0.2647610, 0.0153330, -0.0023320}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
        );
        shells.push_back(
            {
              {14710.0000000, 2207.0000000, 502.8000000, 142.6000000, 46.4700000, 16.7000000, 6.3560000, 1.3160000, 0.3897000},
              {
                {0, false, {-0.0001650, -0.0013080, -0.0064950, -0.0266910, -0.0736900, -0.1707760, -0.1123270, 0.5628140, 0.5687780}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.3897000},
              {
                {0, false, {1.0000000}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {22.6700000, 4.9770000, 1.3470000, 0.3471000},
              {
                {1, false, {0.0448780, 0.2357180, 0.5085210, 0.4581200}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.3471000},
              {
                {1, false, {1.0000000}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.6400000},
              {
                {2, true, {1.0000000}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );

          break;

        default:
          throw "do not know cc-pVDZ basis for this Z";
      }
    }

    return shells;

} */

/*size_t nbasis(const std::vector<libint2::Shell>& shells) {
    size_t n = 0;
    for (const auto& shell: shells)
        n += shell.size();
    return n;
}*/
size_t nbasis(const libint2::BasisSet& basissets) {
  size_t n = 0;
  for (const auto& shell: basissets)
    n += shell.size();
  return n;
}

/*
size_t max_nprim(const std::vector<libint2::Shell>& shells) {
  size_t n = 0;
  for (auto shell: shells)
    n = std::max(shell.nprim(), n);
  return n;
}*/ //### replace by max_nprim()
size_t max_nprim(const std::vector<libint2::BasisSet>& basissets) {
  size_t n = 0;
  for (auto basisset: basissets)
    n = std::max(basisset.max_nprim(), n);
  return n;
}

/*
int max_l(const std::vector<libint2::Shell>& shells) {
  int l = 0;
  for (auto shell: shells)
    for (auto c: shell.contr)
      l = std::max(c.l, l);
  return l;
}*/ //### replace by max_l()
int max_l(const std::vector<libint2::BasisSet>& basissets) {
  int l = 0;
  for (auto basisset: basissets)
      l = std::max(l, int(basisset.max_l()));
    //for (auto c: shell.contr)
      //l = std::max(c.l, l);
  return l;
}

//### Not necessary, replaced by member function
/*
std::vector<size_t> map_shell_to_basis_function(const libint2::BasisSet& basisset) {
  std::vector<size_t> result;
  result.reserve(basisset.size());

  size_t n = 0;
  for (auto basis: basisset) {
    result.push_back(n);
    n += basis.size();
  }

  return result;
}*/

// computes Superposition-Of-Atomic-Densities guess for the molecular density matrix
// in minimal basis; occupies subshells by smearing electrons evenly over the orbitals
Matrix compute_soad(const std::vector<libint2::Atom>& atoms) {

  // compute number of atomic orbitals
  size_t nao = 0;
  for(const auto& atom: atoms) {
    const auto Z = atom.atomic_number;
    if (Z == 1 || Z == 2) // H, He
      nao += 1;
    else if (Z <= 10) // Li - Ne
      nao += 5;
    else
      throw "SOAD with Z > 10 is not yet supported";
  }

  // compute the minimal basis density
  Matrix D = Matrix::Zero(nao, nao);
  size_t ao_offset = 0; // first AO of this atom
  for(const auto& atom: atoms) {
    const auto Z = atom.atomic_number;
    if (Z == 1 || Z == 2) { // H, He
      D(ao_offset, ao_offset) = Z; // all electrons go to the 1s
      ao_offset += 1;
    }
    else if (Z <= 10) {
      D(ao_offset, ao_offset) = 2; // 2 electrons go to the 1s
      D(ao_offset+1, ao_offset+1) = (Z == 3) ? 1 : 2; // Li? only 1 electron in 2s, else 2 electrons
      // smear the remaining electrons in 2p orbitals
      const double num_electrons_per_2p = (Z > 4) ? (double)(Z - 4)/3 : 0;
      for(auto xyz=0; xyz!=3; ++xyz)
        D(ao_offset+2+xyz, ao_offset+2+xyz) = num_electrons_per_2p;
      ao_offset += 5;
    }
  }

  return D * 0.5; // we use densities normalized to # of electrons/2
}

Matrix compute_1body_ints(const libint2::BasisSet& obs,
                          libint2::Operator obtype,
                          const std::vector<libint2::Atom>& atoms)
{
  using libint2::BasisSet;
  //using libint2::Shell;
  using libint2::Engine;
  using libint2::Operator;

  const auto n = nbasis(obs);
  //const auto n = nbasis(shells);
  Matrix result(n,n);

  // construct the overlap integrals engine
  Engine engine(obtype, max_nprim(obs), max_l(obs), 0);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical charges
  if (obtype == Operator::nuclear) {
    std::vector<std::pair<double,std::array<double,3>>> q;
    for(const auto& atom : atoms) {
      q.push_back( {static_cast<double>(atom.atomic_number), {{atom.x, atom.y, atom.z}}} );
    }
    engine.set_params(q);
  }

  //### Not necessary? //commented def
  //auto shell2bf = map_shell_to_basis_function(shells);
  auto shell2bf = obs.shell2bf();

  // buf[0] points to the target shell set after every call  to engine.compute()
  const auto& buf = engine.results();

  // loop over unique shell pairs, {s1,s2} such that s1 >= s2
  // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
  for(auto s1=0; s1!=obs.size(); ++s1) {

	//### not necessary? //###
    auto bf1 = shell2bf[s1]; // first basis function in this shell
    auto n1 = obs[s1].size();

    for(auto s2=0; s2<=s1; ++s2) {

      auto bf2 = shell2bf[s2]; //### Remove?
      auto n2 = obs[s2].size();

      // compute shell pair; return is the pointer to the buffer
      engine.compute(obs[s1], obs[s2]);

      // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
      Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
      result.block(bf2, bf1, n2, n1) = buf_mat.transpose();

    }
  }

  return result;
}

Matrix compute_2body_fock_simple(const libint2::BasisSet& obs,
                                 const Matrix& D) {

  using libint2::BasisSet;
  using libint2::Engine;
  using libint2::Operator;

  const auto n = nbasis(obs);
  //const auto n = nbasis(shells);
  Matrix G = Matrix::Zero(n,n);

  // construct the electron repulsion integrals engine
  Engine engine(Operator::coulomb, /*max_nprim(shells)*/ max_nprim(obs), /*max_l(shells)*/ max_nprim(obs), 0);

  //### not necessary? Commented definintion
  //auto shell2bf = map_shell_to_basis_function(obs);
  auto shell2bf = obs.shell2bf();

  // buf[0] points to the target shell set after every call  to engine.compute()
  const auto& buf = engine.results();

  // loop over shell pairs of the Fock matrix, {s1,s2}
  // Fock matrix is symmetric, but skipping it here for simplicity (see compute_2body_fock)
  for(auto s1=0; s1!=obs.size(); ++s1) {

    auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = obs[s1].size();

    for(auto s2=0; s2!=obs.size(); ++s2) {

      auto bf2_first = shell2bf[s2];
      auto n2 = obs[s2].size();

      // loop over shell pairs of the density matrix, {s3,s4}
      // again symmetry is not used for simplicity
      for(auto s3=0; s3!=obs.size(); ++s3) {

        auto bf3_first = shell2bf[s3];
        auto n3 = obs[s3].size();

        for(auto s4=0; s4!=obs.size(); ++s4) {

          auto bf4_first = shell2bf[s4];
          auto n4 = obs[s4].size();

          // Coulomb contribution to the Fock matrix is from {s1,s2,s3,s4} integrals
          engine.compute(obs[s1], obs[s2], obs[s3], obs[s4]);
          const auto* buf_1234 = buf[0];
          if (buf_1234 == nullptr)
            continue; // if all integrals screened out, skip to next quartet

          // we don't have an analog of Eigen for tensors (yet ... see github.com/BTAS/BTAS, under development)
          // hence some manual labor here:
          // 1) loop over every integral in the shell set (= nested loops over basis functions in each shell)
          // and 2) add contribution from each integral
          for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;
                  G(bf1,bf2) += D(bf3,bf4) * 2.0 * buf_1234[f1234];
                }
              }
            }
          }

          // exchange contribution to the Fock matrix is from {s1,s3,s2,s4} integrals
          engine.compute(obs[s1], obs[s3], obs[s2], obs[s4]);
          const auto* buf_1324 = buf[0];

          for(auto f1=0, f1324=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f3=0; f3!=n3; ++f3) {
              const auto bf3 = f3 + bf3_first;
              for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1324) {
                  const auto bf4 = f4 + bf4_first;
                  G(bf1,bf2) -= D(bf3,bf4) * buf_1324[f1324];
                }
              }
            }
          }

        }
      }
    }
  }

  return G;
}

Matrix compute_2body_fock(const libint2::BasisSet& obs,
                          const Matrix& D) {

  using libint2::BasisSet;
  using libint2::Engine;
  using libint2::Operator;

  std::chrono::duration<double> time_elapsed = std::chrono::duration<double>::zero();

  const auto n = nbasis(obs);
  //const auto n = nbasis(shells);
  Matrix G = Matrix::Zero(n,n);

  // construct the 2-electron repulsion integrals engine
  Engine engine(Operator::coulomb, max_nprim(obs) /*max_nprim(shells)*/, /*max_l(shells)*/ max_l(obs), 0);

  //### Not necessary? commented def
  //auto shell2bf = map_shell_to_basis_function(shells);
  auto shell2bf = obs.shell2bf();

  const auto& buf = engine.results();

  // The problem with the simple Fock builder is that permutational symmetries of the Fock,
  // density, and two-electron integrals are not taken into account to reduce the cost.
  // To make the simple Fock builder efficient we must rearrange our computation.
  // The most expensive step in Fock matrix construction is the evaluation of 2-e integrals;
  // hence we must minimize the number of computed integrals by taking advantage of their permutational
  // symmetry. Due to the multiplicative and Hermitian nature of the Coulomb kernel (and realness
  // of the Gaussians) the permutational symmetry of the 2-e ints is given by the following relations:
  //
  // (12|34) = (21|34) = (12|43) = (21|43) = (34|12) = (43|12) = (34|21) = (43|21)
  //
  // (here we use chemists' notation for the integrals, i.e in (ab|cd) a and b correspond to
  // electron 1, and c and d -- to electron 2).
  //
  // It is easy to verify that the following set of nested loops produces a permutationally-unique
  // set of integrals:
  // foreach a = 0 .. n-1
  //   foreach b = 0 .. a
  //     foreach c = 0 .. a
  //       foreach d = 0 .. (a == c ? b : c)
  //         compute (ab|cd)
  //
  // The only complication is that we must compute integrals over shells. But it's not that complicated ...
  //
  // The real trick is figuring out to which matrix elements of the Fock matrix each permutationally-unique
  // (ab|cd) contributes. STOP READING and try to figure it out yourself. (to check your answer see below)

  // loop over permutationally-unique set of shells
  for(auto s1=0; s1!=obs.size(); ++s1) {

    auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = obs[s1].size();   // number of basis functions in this shell

    for(auto s2=0; s2<=s1; ++s2) {

      auto bf2_first = shell2bf[s2];
      auto n2 = obs[s2].size();

      for(auto s3=0; s3<=s1; ++s3) {

        auto bf3_first = shell2bf[s3];
        auto n3 = obs[s3].size();

        const auto s4_max = (s1 == s3) ? s2 : s3;
        for(auto s4=0; s4<=s4_max; ++s4) {

          auto bf4_first = shell2bf[s4];
          auto n4 = obs[s4].size();

          // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
          auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
          auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
          auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
          auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

          const auto tstart = std::chrono::high_resolution_clock::now();

          engine.compute(obs[s1], obs[s2], obs[s3], obs[s4]);
          const auto* buf_1234 = buf[0];
          if (buf_1234 == nullptr)
            continue; // if all integrals screened out, skip to next quartet

          const auto tstop = std::chrono::high_resolution_clock::now();
          time_elapsed += tstop - tstart;

          // ANSWER
          // 1) each shell set of integrals contributes up to 6 shell sets of the Fock matrix:
          //    F(a,b) += (ab|cd) * D(c,d)
          //    F(c,d) += (ab|cd) * D(a,b)
          //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
          //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
          //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
          //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
          // 2) each permutationally-unique integral (shell set) must be scaled by its degeneracy,
          //    i.e. the number of the integrals/sets equivalent to it
          // 3) the end result must be symmetrized
          for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;

                  const auto value = buf_1234[f1234];

                  const auto value_scal_by_deg = value * s1234_deg;

                  G(bf1,bf2) += D(bf3,bf4) * value_scal_by_deg;
                  G(bf3,bf4) += D(bf1,bf2) * value_scal_by_deg;
                  G(bf1,bf3) -= 0.25 * D(bf2,bf4) * value_scal_by_deg;
                  G(bf2,bf4) -= 0.25 * D(bf1,bf3) * value_scal_by_deg;
                  G(bf1,bf4) -= 0.25 * D(bf2,bf3) * value_scal_by_deg;
                  G(bf2,bf3) -= 0.25 * D(bf1,bf4) * value_scal_by_deg;
                }
              }
            }
          }

        }
      }
    }
  }

  // symmetrize the result and return
  Matrix Gt = G.transpose();
  return 0.5 * (G + Gt);
}
