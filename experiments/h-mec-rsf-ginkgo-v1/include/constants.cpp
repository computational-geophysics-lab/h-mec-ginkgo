#include <ginkgo/ginkgo.hpp>
#include "constants.hpp"





const std::unique_ptr<gko::matrix::Dense<>> x_gko = GinkgoLinSpaced(Nx, xbeg, xend); // (Nx x 1) Vector, horizontal coordinates of basic grid points
const double* x = x_gko->get_values(); // (Nx x 1) Vector, horizontal coordinates of basic grid points
const std::unique_ptr<gko::matrix::Dense<>> y_gko = GinkgoLinSpaced(Ny, ybeg, yend); // (Ny x 1) Vector, vertical coordinates of basic grid points
const double* y = y_gko->get_values(); // (Ny x 1) Vector, vertical coordinates of basic grid points
const std::unique_ptr<gko::matrix::Dense<>> xvx_gko = GinkgoLinSpaced(Nx1, xbeg, xend + dx); // (Nx1 x 1) Vector, Horizontal coordinates of Vx - nodes
const double* xvx = xvx_gko->get_values(); // (Nx1 x 1) Vector, Horizontal coordinates of Vx - nodes
const std::unique_ptr<gko::matrix::Dense<>> yvx_gko = GinkgoLinSpaced(Ny1, ybeg - dy / 2., yend + dy / 2.); // (Ny1 x 1) Vector, Vertical coordinates of Vx - nodes
const double* yvx = yvx_gko->get_values(); // (Ny1 x 1) Vector, Vertical coordinates of Vx - nodes
const std::unique_ptr<gko::matrix::Dense<>> xvy_gko = GinkgoLinSpaced(Nx1, xbeg - dx / 2., xend + dx / 2.); // (Nx1 x 1) Vector, Horizontal coordinates of Vy - nodes
const double* xvy = xvy_gko->get_values(); // (Nx1 x 1) Vector, Horizontal coordinates of Vy - nodes
const std::unique_ptr<gko::matrix::Dense<>> yvy_gko = GinkgoLinSpaced(Ny1, ybeg, yend + dy);// (Ny1 x 1) Vector, Vertical coordinates of Vy - nodes
const double* yvy = yvy_gko->get_values(); // (Ny1 x 1) Vector, Vertical coordinates of Vy - nodes
const std::unique_ptr<gko::matrix::Dense<>> xp_gko = GinkgoLinSpaced(Nx1, xbeg-dx/2., xend+dx/2.); // (Nx x 1) Vector, Horizontal coordinates of P - nodes
const double* xp = xp_gko->get_values(); // (Nx x 1) Vector, Horizontal coordinates of P - nodes
const std::unique_ptr<gko::matrix::Dense<>> yp_gko = GinkgoLinSpaced(Ny1, ybeg-dy/2., yend+dy/2.); // (Ny1 x 1) Vector, Vertical coordinates of P - nodes
const double* yp = yp_gko->get_values(); // (Ny1 x 1) Vector, Vertical coordinates of P - nodes


