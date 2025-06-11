#include "stiffness.h"

RodStiffness computeRodStiffness(const GeomParams& geom, const Material& material) {
    double EA = material.youngs_rod * (geom.axs != nullptr ? *geom.axs : (M_PI * geom.rod_r0 * geom.rod_r0));
    double EI1, EI2;
    if (geom.ixs1 && geom.ixs2) {
        EI1 = material.youngs_rod * (*geom.ixs1);
        EI2 = material.youngs_rod * (*geom.ixs2);
    } else {
        EI1 = EI2 = material.youngs_rod * M_PI * std::pow(geom.rod_r0, 4) / 4.0;
    }
    double GJ = (material.youngs_rod / (2.0 * (1.0 + material.poisson_rod))) *
                (geom.jxs ? *geom.jxs : M_PI * std::pow(geom.rod_r0, 4) / 2.0);

    return {EA, EI1, EI2, GJ};
}

ShellStiffness computeShellStiffness(const GeomParams& geom, const Material& material,
                                     const std::vector<double>& ref_len, bool use_mid_edge) {
    std::vector<double> ks(ref_len.size());
    double kb;

    if (use_mid_edge) {
        double factor = 2.0 * material.youngs_shell * geom.shell_h /
                        (1.0 - std::pow(material.poisson_shell, 2));
        for (size_t i = 0; i < ref_len.size(); ++i) {
            ks[i] = factor * ref_len[i];
        }
        kb = material.youngs_shell * std::pow(geom.shell_h, 3) /
             (24.0 * (1.0 - std::pow(material.poisson_shell, 2)));
    } else {
        double factor = (std::sqrt(3.0) / 2.0) * material.youngs_shell * geom.shell_h;
        for (size_t i = 0; i < ref_len.size(); ++i) {
            ks[i] = factor * ref_len[i];
        }
        kb = (2.0 / std::sqrt(3.0)) * material.youngs_shell * std::pow(geom.shell_h, 3) / 12.0;
    }

    return {ks, kb};
}