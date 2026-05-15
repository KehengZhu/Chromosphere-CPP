#include "../chromosphere.hpp"
#include "armadillo"

namespace chromosphere {

//--------------------------------------------------------------
// Permutation functions

// By Keheng: returns xn(i+1, j, k), a left shift in the x-direction.
// The boundary condition is likely dU/dx = 0 so the last column is copied.

// if xn is a cube, nk = num_of_eq
// if xn is a slice of a cube, nk = 1
Vec ip1(const Vec& xn, const uword nk) {
	// instantiate
	Vec xn_ip1(arma::size(xn));
	// permutate
#pragma omp parallel for collapse(2)
	for(uword i = 0; i < ns - 1; i++) {
		for(uword k = 0; k < nk; k++) {
			xn_ip1(sub2ind(arma::size(ns, num_of_eq), i, k)) =
				xn(sub2ind(arma::size(ns, num_of_eq), i + 1, k));
		}
	}

	// default B.C. (copy the last column, so the BC is likely dU/dx = 0)
	for(uword k = 0; k < nk; k++) {
		xn_ip1(sub2ind(arma::size(ns, num_of_eq), ns - 1, k)) =
			outer_boundary0_i[k];
		if(USE_NEUMANN_BC || nk==SLICE) xn_ip1(sub2ind(arma::size(ns, num_of_eq), ns - 1, k)) = xn(sub2ind(arma::size(ns, num_of_eq), ns-1, k));
	}

	return xn_ip1;
}

// By Keheng: returns xn(i-1, j, k), a right shift in the x-direction.
// The boundary condition is likely dU/dx = 0 so the first column is copied.
Vec im1(const Vec& xn, const uword nk) {
	// instantiate
	Vec xn_im1(arma::size(xn));
	// permutate
#pragma omp parallel for collapse(2)
	for(uword i = 1; i < ns; i++) {
		for(uword k = 0; k < nk; k++) {

			xn_im1(sub2ind(arma::size(ns, num_of_eq), i, k)) =
				xn(sub2ind(arma::size(ns, num_of_eq), i - 1, k));
		}
	}

	// default B.C.
	for(uword k = 0; k < nk; k++) {
		xn_im1(sub2ind(arma::size(ns, num_of_eq), 0, k)) =
			inner_boundary0_i[k];
		if(USE_NEUMANN_BC || nk==SLICE) xn_im1(sub2ind(arma::size(ns, num_of_eq), 0, k)) = xn(sub2ind(arma::size(ns, num_of_eq), 0, k));
	}

	return xn_im1;
}

Vec ip2(const Vec& xn, const uword nk) {
	// instantiate
	Vec xn_ip2(arma::size(xn));
	xn_ip2 = ip1(ip1(xn));

	if(USE_NEUMANN_BC || nk==SLICE) return xn_ip2;

	// default B.C. (copy the last column, so the BC is likely dU/dx = 0)
	for(uword k = 0; k < nk; k++) {
		xn_ip2(sub2ind(arma::size(ns, num_of_eq), ns - 1, k)) =
			outer_boundary1_i[k];
	}

	return xn_ip2;
}

Vec im2(const Vec& xn, const uword nk) {
	// instantiate
	Vec xn_im2(arma::size(xn));
	xn_im2 = im1(im1(xn));

	if(USE_NEUMANN_BC || nk==SLICE) return xn_im2;

	// default B.C.
	for(uword k = 0; k < nk; k++) {
		xn_im2(sub2ind(arma::size(ns, num_of_eq), 0, k)) =
			inner_boundary1_i[k];
	}

	return xn_im2;
}

//--------------------------------------------------------------
// Scalar functions

/* By Keheng:
This is an inversion of scalar_to(). It takes a 3D vector (xn) and extracts the
data from the third dimension (index) and returns it as a 2D scalar field.
*/
Vec get_scalar(const Vec& xn_ii, const uword& index) {
	Vec scalar(ns);

	for(uword i = 0; i < ns; i++) {
		scalar(i) = xn_ii(sub2ind(arma::size(ns, num_of_eq), i, index));
	}

	return scalar;
}

/* By Keheng:
The scalar_to function effectively converts a 2D scalar field (scalar) into a 3D
vector (xn) where the data from scalar is copied along the third dimension
(index). For the rest of the dimensions, the data is zeroed out.

a possibly better name: expand_scalar_to_3D
*/
Vec scalar_to(const Vec& xn_i, const uword& index) {
	Vec xn_ii = zeros<Vec>(ns * num_of_eq);

	for(uword i = 0; i < ns; i++) {
		xn_ii(sub2ind(arma::size(ns, num_of_eq), i, index)) = xn_i(i);
	}

	return xn_ii;
}

//--------------------------------------------------------------
// Operator functions

// Use divided differences to calculate the derivative.
Vec deriv_s(const Vec& xn, const uword type) {
	// forward difference, backward difference, or central difference
	if(type == FORWARD)
		return (ip1(xn, SLICE) - xn) / ds_i;
	else if(type == BACKWARD)
		return (xn - im1(xn, SLICE)) / ds_i;
	return (ip1(xn, SLICE) - im1(xn, SLICE)) / (2.0 * ds_i);
}

// norm of a vector
Vec vec_abs(const Vec& Ax, const Vec& Ay, const Vec& Az) {
	return arma::sqrt(Ax % Ax + Ay % Ay + Az % Az);
}

Vec vec_abs2(const Vec& Ax, const Vec& Ay, const Vec& Az) {
	return Ax % Ax + Ay % Ay + Az % Az;
}

// Convert primitive variables to conserved variables
Vec cons2prim(const Vec& cons) {
	Vec prim(arma::size(cons), fill::zeros);
	const Vec rhoi = get_scalar(cons, CNI);
    const Vec rhon = get_scalar(cons, CNN);
    const Vec rhov = get_scalar(cons, CNV);
    const Vec rhou = get_scalar(cons, CNU);
    const Vec ei = get_scalar(cons, CEI);
    const Vec en = get_scalar(cons, CEN);

    const Vec ni = rhoi/m_i;
    const Vec nn = rhon/m_n;

    const Vec vv = rhov / rhoi;
    const Vec uu = rhou / rhon;
    const Vec phig = 0.5*(gPotential_imh+gPotential_iph);
    const Vec pi = 2.0/3.0*ei - 1.0/3.0*rhoi%vv%vv - 2.0/3.0*rhoi%phig; 
    const Vec pn = 2.0/3.0*en - 1.0/3.0*rhon%uu%uu - 2.0/3.0*rhon%phig;

	prim += scalar_to(rhoi, PNI);
	prim += scalar_to(rhon, PNN);
	prim += scalar_to(vv, PV);
	prim += scalar_to(uu, PU);
	prim += scalar_to(pi, PPI);
	prim += scalar_to(pn, PPN);
	return prim;
}

Vec prim2cons(const Vec& prim) {
	Vec cons(arma::size(prim), fill::zeros);
	const Vec rhoi = get_scalar(prim, PNI);
	const Vec rhon = get_scalar(prim, PNN);
	const Vec vv = get_scalar(prim, PV);
	const Vec uu = get_scalar(prim, PU);
	const Vec pi = get_scalar(prim, PPI);
	const Vec pn = get_scalar(prim, PPN);

	const Vec rhov = rhoi % vv;
	const Vec rhou = rhon % uu;
	const Vec phig = 0.5*(gPotential_imh+gPotential_iph);
	const Vec ei = 3.0/2.0*pi + 0.5*rhoi%vv%vv + rhoi%phig;
	const Vec en = 3.0/2.0*pn + 0.5*rhon%uu%uu + rhon%phig;

	cons += scalar_to(rhoi, CNI);
	cons += scalar_to(rhon, CNN);
	cons += scalar_to(rhov, CNV);
	cons += scalar_to(rhou, CNU);
	cons += scalar_to(ei, CEI);
	cons += scalar_to(en, CEN);
	return cons;
}

Vec shock_tube_ic() {
	Vec p_xn = zeros<Vec>(ns * num_of_eq);
	// const float b = 0.75, b_l = 1.0, n_l = 1.0, p_l = 1.0, b_r = -1.0, n_r = 0.125, p_r = 0.1;
	const float b = 0, b_l = 0, n_l = 1.0, p_l = 1.0, b_r = 0, n_r = 0.125, p_r = 0.1;
#pragma omp parallel for collapse(2)
	for(uword i = 0; i < ns; i++) {
		uword ij_ni = sub2ind(arma::size(ns, num_of_eq), i, PNI);
		uword ij_nn = sub2ind(arma::size(ns, num_of_eq), i, PNN);
		uword ij_Pi = sub2ind(arma::size(ns, num_of_eq), i, PPI);
		uword ij_Pn = sub2ind(arma::size(ns, num_of_eq), i, PPN);
		if(i > ns / 2) {
			p_xn(ij_ni) = n_r;
			p_xn(ij_nn) = n_r;
			p_xn(ij_Pi) = p_r;
			p_xn(ij_Pn) = p_r;
		}
		else {
			p_xn(ij_ni) = n_l;
			p_xn(ij_nn) = n_l;
			p_xn(ij_Pi) = p_l;
			p_xn(ij_Pn) = p_l;
		}
	}
	return prim2cons(p_xn);
}

}; // end namespace chromosphere

void array2vec(const float* from_array, Vec& to_vec, const uword& size) {
	for(uword i = 0; i < size; i++) {
		to_vec(i) = from_array[i];
	}
}

void vec2array(const Vec& from_vec, float* to_array, const uword& size) {
	for(uword i = 0; i < size; i++) {
		to_array[i] = from_vec(i);
	}
}
