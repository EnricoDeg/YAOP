/* Copyright (C) 2023  Enrico Degregori, Wilton Jaciel Loch
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <memory>
#include <cstdlib>
#include <fstream>

#include "src/TKE.hpp"

template<typename T>
void read_input(T *data, const std::string &filename, int npoints, int nproma, int nlevs, int nblocks, int npromz);

template<typename T>
void init_fix(T *data, int nproma, int nlevs, int nblocks, T value);

template<typename T>
void read_input_raw(T *data, const std::string &filename, int npoints);

template<typename T>
void write_output(T*data, const std::string &filename, int npoints);

int main(int argc, char ** argv) {
    int nproma = 15117;
    int nlevs = 40;
    int ntimesteps = 10;
    int ncells = 15105;
    int nedges = 23207;

    int nblocks_cells = ncells / nproma + 1;
    int npromz_cells  = ncells % nproma;
    int nblocks_edges = nedges / nproma + 1;
    int npromz_edges  = nedges % nproma;

    int edges_block_size = nproma;
    int edges_start_block = 0;
    int edges_end_block = nblocks_edges - 1;
    int edges_start_index = 0;
    int edges_end_index = npromz_edges - 1;
    int cells_block_size = nproma;
    int cells_start_block = 0;
    int cells_end_block = nblocks_cells - 1;
    int cells_start_index = 0;
    int cells_end_index = npromz_cells - 1;

    int vert_mix_type = 2;
    int vmix_idemix_tke = 4;
    int vert_cor_type = 0;
    double dtime = 0.0;
    double OceanReferenceDensity = 1025.022;
    double grav = 9.80665;
    int l_lc = 0;
    double clc = 0.15;
    double ReferencePressureIndbars = 1035.0*grav*1.0e-4;
    double pi = 3.14159265358979323846264338327950288;

    std::shared_ptr<TKE> ocean_physics;
    // Initialize TKE
    ocean_physics.reset(new TKE(nproma, nlevs, nblocks_cells, vert_mix_type, vmix_idemix_tke,
                                vert_cor_type, dtime, OceanReferenceDensity, grav,
                                l_lc, clc, ReferencePressureIndbars, pi));

    // Allocate arrays
    double *depth_CellInterface = reinterpret_cast<double *>(
                                  malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(depth_CellInterface, "examples/input/depth_CellInterface",
                       ncells, nproma, nlevs+1, nblocks_cells, npromz_cells);

    double *prism_center_dist_c = reinterpret_cast<double *>(
                                  malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(prism_center_dist_c, "examples/input/prism_center_dist_c",
                       ncells, nproma, nlevs+1, nblocks_cells, npromz_cells);

    double *inv_prism_center_dist_c = reinterpret_cast<double *>(
                                      malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(inv_prism_center_dist_c, "examples/input/inv_prism_center_dist_c",
                       ncells, nproma, nlevs+1, nblocks_cells, npromz_cells);

    double *prism_thick_c = reinterpret_cast<double *>(malloc(nproma * nlevs * nblocks_cells * sizeof(double)));
    read_input<double>(prism_thick_c, "examples/input/prism_thick_c",
                       ncells, nproma, nlevs, nblocks_cells, npromz_cells);

    int *dolic_c = reinterpret_cast<int *>(malloc(nproma * nblocks_cells * sizeof(int)));
    read_input<int>(dolic_c, "examples/input/dolic_c", ncells, nproma, 1, nblocks_cells, npromz_cells);

    int *dolic_e = reinterpret_cast<int *>(malloc(nproma * nblocks_edges * sizeof(int)));
    read_input_raw<int>(dolic_e, "dolic_e", nproma*nblocks_edges);

    double *zlev_i = reinterpret_cast<double *>(malloc(nlevs * sizeof(double)));
    read_input<double>(zlev_i, "examples/input/zlev_i", 1, 1, nlevs, 1, 1);

    double *wet_c = reinterpret_cast<double *>(malloc(nproma * nlevs * nblocks_cells * sizeof(double)));
    read_input<double>(wet_c, "examples/input/wet_c", ncells, nproma, nlevs, nblocks_cells, npromz_cells);

    int *edges_cell_idx = reinterpret_cast<int *>(malloc(nproma * 2 * nblocks_edges * sizeof(int)));
    read_input_raw<int>(edges_cell_idx, "edges_cell_idx", 2*nproma*nblocks_edges);

    int *edges_cell_blk = reinterpret_cast<int *>(malloc(nproma * 2 * nblocks_edges * sizeof(int)));
    read_input_raw<int>(edges_cell_blk, "edges_cell_blk", 2*nproma*nblocks_edges);

    double *temp = reinterpret_cast<double *>(malloc(nproma * nlevs * nblocks_cells * sizeof(double)));
    read_input<double>(temp, "examples/input/to", ncells, nproma, nlevs, nblocks_cells, npromz_cells);

    double *salt = reinterpret_cast<double *>(malloc(nproma * nlevs * nblocks_cells * sizeof(double)));
    read_input<double>(salt, "examples/input/so", ncells, nproma, nlevs, nblocks_cells, npromz_cells);

    double *stretch_c = reinterpret_cast<double *>(malloc(nproma * nblocks_cells * sizeof(double)));
    init_fix<double>(stretch_c, nproma, 1, nblocks_cells, 1.0);

    double *eta_c = reinterpret_cast<double *>(malloc(nproma * nblocks_cells * sizeof(double)));
    init_fix<double>(eta_c, nproma, 1, nblocks_cells, 0.0);

    double *p_vn_x1 = reinterpret_cast<double *>(malloc(nproma * nlevs * nblocks_cells * sizeof(double)));
    init_fix<double>(p_vn_x1, nproma, nlevs, nblocks_cells, 0.0);

    double *p_vn_x2 = reinterpret_cast<double *>(malloc(nproma * nlevs * nblocks_cells * sizeof(double)));
    init_fix<double>(p_vn_x2, nproma, nlevs, nblocks_cells, 0.0);

    double *p_vn_x3 = reinterpret_cast<double *>(malloc(nproma * nlevs * nblocks_cells * sizeof(double)));
    init_fix<double>(p_vn_x3, nproma, nlevs, nblocks_cells, 0.0);

    double *tke = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(tke, "examples/input/tke", ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *tke_plc_in = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));

    double *hlc_in = reinterpret_cast<double *>(malloc(nproma * nblocks_cells * sizeof(double)));

    double *wlc_in = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));

    double *u_stokes_in = reinterpret_cast<double *>(malloc(nproma * nblocks_cells * sizeof(double)));

    double *a_veloc_v = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_edges * sizeof(double)));

    double *a_temp_v = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(a_temp_v, "examples/input/A_tracer_v_to",
                       ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *a_salt_v = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(a_salt_v, "examples/input/A_tracer_v_so",
                       ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *iwe_Tdis = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));

    double *cvmix_dummy_1 = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(cvmix_dummy_1, "examples/input/cvmix_dummy_1",
                       ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *cvmix_dummy_2 = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(cvmix_dummy_2, "examples/input/cvmix_dummy_2",
                       ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *cvmix_dummy_3 = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(cvmix_dummy_3, "examples/input/cvmix_dummy_3",
                       ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *tke_Tbpr = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(tke_Tbpr, "examples/input/tke_Tbpr", ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *tke_Tspr = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(tke_Tspr, "examples/input/tke_Tspr", ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *tke_Tdif = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(tke_Tdif, "examples/input/tke_Tdif", ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *tke_Tdis = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(tke_Tdis, "examples/input/tke_Tdis", ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *tke_Twin = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(tke_Twin, "examples/input/tke_Twin", ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *tke_Tiwf = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(tke_Tiwf, "examples/input/tke_Tiwf", ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *tke_Tbck = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(tke_Tbck, "examples/input/tke_Tbck", ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *tke_Ttot = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(tke_Ttot, "examples/input/tke_Ttot", ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *tke_Lmix = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(tke_Lmix, "examples/input/tke_Lmix", ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *tke_Pr = reinterpret_cast<double *>(malloc(nproma * (nlevs+1) * nblocks_cells * sizeof(double)));
    read_input<double>(tke_Pr, "examples/input/tke_Pr", ncells, nproma, (nlevs+1), nblocks_cells, npromz_cells);

    double *stress_xw = reinterpret_cast<double *>(malloc(nproma * nblocks_cells * sizeof(double)));
    read_input<double>(stress_xw, "examples/input/atmos_fluxes_stress_xw",
                       ncells, nproma, 1, nblocks_cells, npromz_cells);

    double *stress_yw = reinterpret_cast<double *>(malloc(nproma * nblocks_cells * sizeof(double)));
    read_input<double>(stress_yw, "examples/input/atmos_fluxes_stress_yw",
                       ncells, nproma, 1, nblocks_cells, npromz_cells);

    double *fu10 = reinterpret_cast<double *>(malloc(nproma * nblocks_cells * sizeof(double)));
    read_input<double>(fu10, "examples/input/Wind_Speed_10m", ncells, nproma, 1, nblocks_cells, npromz_cells);

    double *concsum = reinterpret_cast<double *>(malloc(nproma * nblocks_cells * sizeof(double)));
    read_input<double>(concsum, "examples/input/conc", ncells, nproma, 1, nblocks_cells, npromz_cells);

    #pragma acc enter data copyin(depth_CellInterface[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(prism_center_dist_c[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(inv_prism_center_dist_c[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(prism_thick_c[0:nproma*nlevs*nblocks_cells-1])
    #pragma acc enter data copyin(dolic_c[0:nproma*nblocks_cells-1], dolic_e[0:nproma*nblocks_edges-1])
    #pragma acc enter data copyin(zlev_i[0:nlevs-1], wet_c[0:nproma*nlevs*nblocks_cells-1])
    #pragma acc enter data copyin(edges_cell_idx[0:nproma*2*nblocks_edges-1])
    #pragma acc enter data copyin(edges_cell_blk[0:nproma*2*nblocks_edges-1])
    #pragma acc enter data copyin(tke[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(tke_plc_in[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(hlc_in[0:nproma*nblocks_cells-1])
    #pragma acc enter data copyin(wlc_in[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(u_stokes_in[0:nproma*nblocks_cells-1])
    #pragma acc enter data copyin(a_veloc_v[0:nproma*(nlevs+1)*nblocks_edges-1])
    #pragma acc enter data copyin(a_temp_v[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(a_salt_v[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(iwe_Tdis[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(cvmix_dummy_1[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(cvmix_dummy_2[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(cvmix_dummy_3[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(tke_Tbpr[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(tke_Tspr[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(tke_Tdif[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(tke_Tdis[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(tke_Twin[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(tke_Tiwf[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(tke_Tbck[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(tke_Ttot[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(tke_Lmix[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(tke_Pr[0:nproma*(nlevs+1)*nblocks_cells-1])
    #pragma acc enter data copyin(temp[0:nproma*nlevs*nblocks_cells-1], salt[0:nproma*nlevs*nblocks_cells-1])
    #pragma acc enter data copyin(stretch_c[0:nproma*nblocks_cells-1], eta_c[0:nproma*nblocks_cells-1])
    #pragma acc enter data copyin(p_vn_x1[0:nproma*nlevs*nblocks_cells-1])
    #pragma acc enter data copyin(p_vn_x2[0:nproma*nlevs*nblocks_cells-1])
    #pragma acc enter data copyin(p_vn_x3[0:nproma*nlevs*nblocks_cells-1])
    #pragma acc enter data copyin(stress_xw[0:nproma*nblocks_cells-1], stress_yw[0:nproma*nblocks_cells-1])
    #pragma acc enter data copyin(fu10[0:nproma*nblocks_cells-1])
    #pragma acc enter data copyin(concsum[0:nproma*nblocks_cells-1])

    for (int t = 0; t < ntimesteps; t++) {
      #pragma acc host_data use_device(depth_CellInterface, prism_center_dist_c, inv_prism_center_dist_c)
      #pragma acc host_data use_device(prism_thick_c, dolic_c, zlev_i, edges_cell_idx, edges_cell_blk, dolic_e, wet_c)
      #pragma acc host_data use_device(tke, tke_plc_in, hlc_in, wlc_in, u_stokes_in, a_veloc_v, a_temp_v, a_salt_v)
      #pragma acc host_data use_device(iwe_Tdis, cvmix_dummy_1, cvmix_dummy_2, cvmix_dummy_3, tke_Tbpr, tke_Tspr)
      #pragma acc host_data use_device(tke_Tdif, tke_Tdis, tke_Twin, tke_Tiwf, tke_Tbck, tke_Ttot, tke_Lmix, tke_Pr)
      #pragma acc host_data use_device(temp, salt, stretch_c, eta_c, stress_xw, stress_yw, fu10, concsum)
      #pragma acc host_data use_device(p_vn_x1, p_vn_x2, p_vn_x3)
      ocean_physics->calc(depth_CellInterface, prism_center_dist_c,
                          inv_prism_center_dist_c, prism_thick_c,
                          dolic_c, dolic_e, zlev_i, wet_c,
                          edges_cell_idx, edges_cell_blk,
                          temp, salt, stretch_c, eta_c,
                          p_vn_x1, p_vn_x2, p_vn_x3,
                          tke, tke_plc_in, hlc_in, wlc_in,
                          u_stokes_in, a_veloc_v, a_temp_v, a_salt_v,
                          iwe_Tdis, cvmix_dummy_1, cvmix_dummy_2,
                          cvmix_dummy_3, tke_Tbpr, tke_Tspr,
                          tke_Tdif, tke_Tdis, tke_Twin,
                          tke_Tiwf, tke_Tbck, tke_Ttot,
                          tke_Lmix, tke_Pr, stress_xw,
                          stress_yw, fu10, concsum,
                          edges_block_size, edges_start_block, edges_end_block,
                          edges_start_index, edges_end_index, cells_block_size,
                          cells_start_block, cells_end_block, cells_start_index,
                          cells_end_index);
      #pragma acc wait
    }

//    #pragma acc update host(tke[0:nproma*(nlevs+1)*nblocks_cells-1])
//    #pragma acc update host(a_temp_v[0:nproma*(nlevs+1)*nblocks_cells-1])
//    #pragma acc update host(a_salt_v[0:nproma*(nlevs+1)*nblocks_cells-1])
//    write_output<double>(tke, "examples/output/tke", nproma*(nlevs+1)*nblocks_cells);
//    write_output<double>(a_temp_v, "examples/output/a_temp_v", nproma*(nlevs+1)*nblocks_cells);
//    write_output<double>(a_salt_v, "examples/output/a_salt_v", nproma*(nlevs+1)*nblocks_cells);

    ocean_physics.reset();

    // Deallocate arrays on host and device
    #pragma acc exit data delete(depth_CellInterface, prism_center_dist_c, inv_prism_center_dist_c)
    #pragma acc exit data delete(prism_thick_c, dolic_c, dolic_e, zlev_i, wet_c, edges_cell_idx, edges_cell_blk)
    #pragma acc exit data delete(tke, tke_plc_in, hlc_in, wlc_in, u_stokes_in, a_veloc_v, a_temp_v, a_salt_v)
    #pragma acc exit data delete(iwe_Tdis, cvmix_dummy_1, cvmix_dummy_2, cvmix_dummy_3, tke_Tbpr, tke_Tspr)
    #pragma acc exit data delete(tke_Tdif, tke_Tdis, tke_Twin, tke_Tiwf, tke_Tbck, tke_Ttot, tke_Lmix, tke_Pr)
    #pragma acc exit data delete(temp, salt, stretch_c, eta_c)
    #pragma acc exit data delete(p_vn_x1, p_vn_x2, p_vn_x3)
    #pragma acc exit data delete(stress_xw, stress_yw)
    #pragma acc exit data delete(fu10)
    #pragma acc exit data delete(concsum)

    free(depth_CellInterface);
    free(prism_center_dist_c);
    free(inv_prism_center_dist_c);
    free(prism_thick_c);
    free(dolic_c);
    free(dolic_e);
    free(zlev_i);
    free(wet_c);
    free(edges_cell_idx);
    free(edges_cell_blk);

    free(tke);
    free(tke_plc_in);
    free(hlc_in);
    free(wlc_in);
    free(u_stokes_in);
    free(a_veloc_v);
    free(a_temp_v);
    free(a_salt_v);
    free(iwe_Tdis);
    free(cvmix_dummy_1);
    free(cvmix_dummy_2);
    free(cvmix_dummy_3);
    free(tke_Tbpr);
    free(tke_Tspr);
    free(tke_Tdif);
    free(tke_Tdis);
    free(tke_Twin);
    free(tke_Tiwf);
    free(tke_Tbck);
    free(tke_Ttot);
    free(tke_Lmix);
    free(tke_Pr);

    free(temp);
    free(salt);
    free(stretch_c);
    free(eta_c);
    free(p_vn_x1);
    free(p_vn_x2);
    free(p_vn_x3);

    free(stress_xw);
    free(stress_yw);

    free(fu10);

    free(concsum);

    return 0;
}

template<typename T>
void read_input(T *data, const std::string &filename, int npoints, int nproma, int nlevs, int nblocks, int npromz) {
    T *buffer = reinterpret_cast<T *>(malloc(npoints * nlevs * sizeof(T)));
    std::ifstream ifile;
    ifile.open(filename);
    int k = 0;
    while (ifile >> buffer[k]) k++;
    ifile.close();

    for (int jb = 0; jb < nblocks; jb++) {
        for (int level = 0; level < nlevs; level++) {
            int nc = nproma;
            if (jb == nblocks-1) nc = npromz;
            for (int jc = 0; jc < nc; jc++) {
                data[jc+level*nproma+jb*nproma*nlevs] = buffer[jc+nproma*jb+level*npoints];
            }
            for (int jc = nc; jc < nproma; jc++)
                data[jc+level*nproma+jb*nproma*nlevs] = 0.0;
        }
    }

    free(buffer);
}

template<typename T>
void init_fix(T *data, int nproma, int nlevs, int nblocks, T value) {
    for (int jb = 0; jb < nblocks; jb++)
        for (int level = 0; level < nlevs; level++)
            for (int jc = 0; jc < nproma; jc++)
                data[jc+level*nproma+jb*nproma*nlevs] = value;
}

template<typename T>
void read_input_raw(T *data, const std::string &filename, int npoints) {
    std::ifstream ifile;
    ifile.open(filename);
    int k = 0;
    while (ifile >> data[k]) k++;
    ifile.close();
    for (int i = k; i < npoints; i++)
        data[i] = 0;
}

template<typename T>
void write_output(T*data, const std::string &filename, int npoints) {
    std::ofstream ofile;
    ofile.open(filename);
    for (int k = 0; k < npoints; k++)
        ofile << data[k] << std::endl;
    ofile.close();
}
