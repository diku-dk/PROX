#include <simulators_prox.h>
#include <simulators_prox_data.h>

#include <util_profiling.h>
#include <util_matlab_write_profiling.h>
#include <util_python_write_profiling.h>
#include <util_python_write_matrix.h>
#include <util_log.h>

#include <prox_rigid_body.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <fstream>

#include <string>
#include <vector>
#include <fstream>   // needed for std::ofstream

namespace simulators
{

bool ProxEngine::write_profilingMatlab(std::string const& filename)
{
    util::Log logging;

    std::string const newline = util::Log::newline();

    std::ofstream matlab;

    matlab.open(filename.c_str(), std::ios::out);

    if (!matlab.is_open())
    {
        logging << "ProxEngine::write_profiling(): error could not open file = " << filename.c_str()
                << util::Log::newline();

        return false;
    }
    return true;
}

bool ProxEngine::write_profiling(std::string const& filename)
{
    util::Log logging;

    std::string const newline = util::Log::newline();

    std::ofstream matlab;

    matlab.open(filename.c_str(), std::ios::out);

    if (!matlab.is_open())
    {
        logging << "ProxEngine::write_profiling(): error could not open file = " << filename.c_str()
                << util::Log::newline();

        return false;
    }

    matlab << "%close all;" << std::endl;
    matlab << "%clear all;" << std::endl;
    matlab << "%clc;" << std::endl;
    matlab << std::endl;
    matlab << std::endl;
    matlab << util::matlab_write_profiling() << std::endl;
    matlab << std::endl;
    matlab << std::endl;

    matlab << std::endl;
    matlab << std::endl;
    matlab << "%red   = [0.7,0.1,0.1];" << std::endl;
    matlab << "%green = [0.1,0.7,0.1];" << std::endl;
    matlab << "%blue  = [0.1,0.1,0.7];" << std::endl;
    matlab << "%f_type = 'Times';" << std::endl;
    matlab << "%f_size = 20;" << std::endl;
    matlab << std::endl;
    matlab << std::endl;
    matlab << "%frames = [1:length(Ekin)];" << std::endl;
    matlab << "%Emech = Epot + Ekin;" << std::endl;
    matlab << std::endl;
    matlab << std::endl;
    matlab << "%figure(100);" << std::endl;
    matlab << "%e1 = plot(frames, Ekin,'r-'); " << std::endl;
    matlab << "%hold on; " << std::endl;
    matlab << "%e2 = plot(frames, Epot, 'g-'); " << std::endl;
    matlab << "%e3 = plot(frames, Emech, 'b-'); " << std::endl;
    matlab << "%hold off; " << std::endl;
    matlab << "%title('Kinetic and potential energy of total system', 'FontSize', f_size, 'FontName', f_type); "
           << std::endl;
    matlab << "%legend([e1,e2,e3], {'Kinetic', 'Potential','Mechanical'}, 'FontSize', f_size, 'FontName', f_type); "
           << std::endl;
    matlab << "%ylabel('Energy (Joules)', 'FontSize', f_size, 'FontName', f_type)" << std::endl;
    matlab << "%xlabel('Frame', 'FontSize', f_size, 'FontName', f_type)" << std::endl;
    matlab << "%print(gcf,'-depsc2','energy');" << std::endl;
    matlab << "%print(gcf,'-dpng','energy');" << std::endl;
    matlab << std::endl;
    matlab << std::endl;
    matlab << "%figure(200);" << std::endl;
    matlab << "%plot(abs(average_penetration)*100);" << std::endl;
    matlab << "%title('Penetrations', 'FontSize', f_size, 'FontName', f_type); " << std::endl;
    matlab << "%xlabel('Frame', 'FontSize', f_size, 'FontName', f_type); " << std::endl;
    matlab << "%ylabel('Penetration depth (cm)', 'FontSize', f_size, 'FontName', f_type)" << std::endl;
    matlab << "%print(gcf,'-depsc2','penetration');" << std::endl;
    matlab << "%print(gcf,'-dpng','penetration');" << std::endl;
    matlab << std::endl;
    matlab << std::endl;
    matlab << "%figure(300);" << std::endl;
    matlab << "%plot(contacts);" << std::endl;
    matlab << "%title('Contacts', 'FontSize', f_size, 'FontName', f_type); " << std::endl;
    matlab << "%xlabel('Frame', 'FontSize', f_size, 'FontName', f_type); " << std::endl;
    matlab << "%ylabel('Number of contacts', 'FontSize', f_size, 'FontName', f_type)" << std::endl;
    matlab << "%print(gcf,'-depsc2','contacts');" << std::endl;
    matlab << "%print(gcf,'-dpng','contacts');" << std::endl;
    matlab << std::endl;
    matlab << std::endl;
    matlab << "%figure(400);" << std::endl;
    matlab << "%plot(contacts(contacts>0) ,solver, '.'); " << std::endl;
    matlab << "%title('Time used by solver vs number of contacts', 'FontSize', f_size, 'FontName', f_type); "
           << std::endl;
    matlab << "%xlabel('Number of contacts', 'FontSize', f_size, 'FontName', f_type); " << std::endl;
    matlab << "%ylabel('Solver (ms)', 'FontSize', f_size, 'FontName', f_type)" << std::endl;
    matlab << "%print(gcf,'-depsc2','solver_time_vs_contacts');" << std::endl;
    matlab << "%print(gcf,'-dpng','solver_time_vs_contacts');" << std::endl;
    matlab << std::endl;
    matlab << std::endl;
    matlab << "%figure(500);" << std::endl;
    matlab << "%plot(frames(contacts>0), solver); " << std::endl;
    matlab << "%title('Time used by solver per frame', 'FontSize', f_size, 'FontName', f_type)" << std::endl;
    matlab << "%xlabel('Frame', 'FontSize', f_size, 'FontName', f_type); " << std::endl;
    matlab << "%ylabel('Solver (ms)', 'FontSize', f_size, 'FontName', f_type)" << std::endl;
    matlab << "%print(gcf,'-depsc2','solver_time_per_frame');" << std::endl;
    matlab << "%print(gcf,'-dpng','solver_time_per_frame');" << std::endl;
    matlab << std::endl;
    matlab << std::endl;
    matlab << "%figure(600);" << std::endl;
    matlab << "%semilogy(cell2mat(convergence(1)));" << std::endl;
    matlab << "%hold on" << std::endl;
    matlab << "%for frame = 2:numel(convergence); " << std::endl;
    matlab << "%  semilogy(cell2mat(convergence(frame)));" << std::endl;
    matlab << "%end" << std::endl;
    matlab << "%title('Convergence', 'FontSize', f_size, 'FontName',f_type );" << std::endl;
    matlab << "%xlabel('Solver iteration', 'FontSize', f_size, 'FontName', f_type); " << std::endl;
    matlab << "%ylabel('Natural merit function', 'FontSize', f_size, 'FontName', f_type);" << std::endl;
    matlab << "%print(gcf,'-depsc2','convergence');" << std::endl;
    matlab << "%print(gcf,'-dpng','convergence');" << std::endl;
    matlab << std::endl;
    matlab << std::endl;
    matlab << "%figure(700);" << std::endl;
    matlab << "%plot(cell2mat(rfactor(1)));" << std::endl;
    matlab << "%hold on" << std::endl;
    matlab << "%for frame = 2:numel(rfactor);" << std::endl;
    matlab << "%  plot(cell2mat(rfactor(frame)));" << std::endl;
    matlab << "%end" << std::endl;
    matlab << "%title('R-factor development', 'FontSize', f_size, 'FontName',f_type );" << std::endl;
    matlab << "%xlabel('Solver iteration', 'FontSize', f_size, 'FontName', f_type); ";
    matlab << "%ylabel('R-factor', 'FontSize', f_size, 'FontName', f_type);";
    matlab << "%print(gcf,'-depsc2','rfactors');" << std::endl;
    matlab << "%print(gcf,'-dpng','rfactors');" << std::endl;
    matlab << std::endl;
    matlab << std::endl;
    matlab.flush();
    matlab.close();

    logging << "ProxEngine::write_profiling(): Done writing profile data..." << newline;

    return true;
}

bool ProxEngine::write_profiling_python(std::string const& filename)
{
    util::Log logging;

    std::string const newline = util::Log::newline();

    std::ofstream python;

    python.open(filename.c_str(), std::ios::out);

    if (!python.is_open())
    {
        logging << "ProxEngine::write_profiling(): error could not open file = " << filename.c_str()
                << util::Log::newline();

        return false;
    }

    python << std::endl;
    python << std::endl;
    python << util::python_write_profiling() << std::endl;
    python << std::endl;
    python << std::endl;

    python << "import numpy as np" << "\n";
    python << "import matplotlib.pyplot as plt" << "\n";
    python << "import matplotlib as mpl" << "\n";
    python << "mpl.rcParams['lines.linewidth'] = 0.75" << "\n";
    python << "mpl.rcParams['lines.markersize'] = 6" << "\n";
    python << "\n";
    python << "f_type = 'Times';" << "\n";
    python << "f_size = 20;" << "\n";
    python << "\n";
    python << "\n";
    python << "import numpy as np" << "\n";
    python << "import matplotlib.pyplot as plt" << "\n";
    python << "import warnings" << "\n";
    python << "\n";
    python << "def toNp(x):" << "\n";
    python << "    return np.asarray(x, dtype=float)" << "\n";
    python << "\n";
    python << "\n";
    python << "Ekin = toNp(Ekin)" << "\n";
    python << "Epot = toNp(Epot)" << "\n";
    python << "Emech = Epot + Ekin" << "\n";
    python << "frames = np.arange(1, Ekin.size + 1)" << "\n";
    python << "\n";
    python << "\n";
    python << "plt.figure(100)" << "\n";
    python << "e1, = plt.plot(frames, Ekin, 'r-', label='Kinetic')" << "\n";
    python << "e2, = plt.plot(frames, Epot, 'g-', label='Potential')" << "\n";
    python << "e3, = plt.plot(frames, Emech, 'b-', label='Mechanical')" << "\n";
    python << "plt.title('Kinetic and potential energy of total system', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.legend(handles=[e1, e2, e3], fontsize=f_size)" << "\n";
    python << "plt.ylabel('Energy (Joules)', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.xlabel('Frame', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.grid(True, linestyle=':', linewidth=0.5)" << "\n";
    python << "plt.tight_layout()" << "\n";
    python << "plt.savefig('energy_py.eps', format='eps', bbox_inches='tight')" << "\n";
    python << "plt.savefig('energy_py.png', dpi=300, bbox_inches='tight')" << "\n";
    python << "\n";
    python << "\n";
    python << "plt.figure(200)" << "\n";
    python << "avg_pen = toNp(average_penetration)" << "\n";
    python << "plt.plot(np.abs(avg_pen) * 100)" << "\n";
    python << "plt.title('Penetrations', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.xlabel('Frame', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.ylabel('Penetration depth (cm)', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.grid(True, linestyle=':', linewidth=0.5)" << "\n";
    python << "plt.tight_layout()" << "\n";
    python << "plt.savefig('penetration_py.eps', format='eps', bbox_inches='tight')" << "\n";
    python << "plt.savefig('penetration_py.png', dpi=300, bbox_inches='tight')" << "\n";
    python << "\n";
    python << "\n";
    python << "plt.figure(300)" << "\n";
    python << "plt.plot(toNp(contacts))" << "\n";
    python << "plt.title('Contacts', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.xlabel('Frame', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.ylabel('Number of contacts', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.grid(True, linestyle=':', linewidth=0.5)" << "\n";
    python << "\n";
    python << "plt.tight_layout()" << "\n";
    python << "plt.savefig('contacts_py.eps', format='eps', bbox_inches='tight')" << "\n";
    python << "plt.savefig('contacts_py.png', dpi=300, bbox_inches='tight')" << "\n";
    python << "\n";
    python << "\n";
    python << "\n";
    python << "plt.figure(400)" << "\n";
    python << "contactsNp = toNp(contacts)" << "\n";
    python << "solverNp = toNp(solver)" << "\n";
    python << "\n";
    python << "mask = contactsNp > 0" << "\n";
    python << "xVals = contactsNp[mask]" << "\n";
    python << "\n";
    python << "if solverNp.size == xVals.size:" << "\n";
    python << "    yVals = solverNp" << "\n";
    python << "elif solverNp.size == contactsNp.size:" << "\n";
    python << "    yVals = solverNp[mask]" << "\n";
    python << "else:" << "\n";
    python << "    minLen = min(solverNp.size, xVals.size)" << "\n";
    python << "    xVals = xVals[:minLen]" << "\n";
    python << "    yVals = solverNp[:minLen]" << "\n";
    python << "\n";
    python << "plt.plot(xVals, yVals, '.', markersize=6)" << "\n";
    python << "plt.title('Time used by solver vs number of contacts', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.xlabel('Number of contacts', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.ylabel('Solver (ms)', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.grid(True, linestyle=':', linewidth=0.5)" << "\n";
    python << "plt.tight_layout()" << "\n";
    python << "\n";
    python << "plt.savefig('solver_time_vs_contacts_py.eps', format='eps', bbox_inches='tight')" << "\n";
    python << "plt.savefig('solver_time_vs_contacts_py.png', dpi=300, bbox_inches='tight')" << "\n";
    python << "\n";
    python << "\n";
    python << "\n";
    python << "plt.figure(500)" << "\n";
    python << "if solverNp.size == frames.size:" << "\n";
    python << "    xVal = frames[mask]" << "\n";
    python << "    yVal = solverNp[mask]" << "\n";
    python << "elif solverNp.size == mask.sum():" << "\n";
    python << "    xVal = frames[mask]" << "\n";
    python << "    yVal = solverNp" << "\n";
    python << "else:" << "\n";
    python << "    minLen = min(frames[mask].size, solverNp.size)" << "\n";
    python << "    print(f'warning: solver length ({solverNp.size}) not matching frames(mask) length "
              "({frames[mask].size}). Truncating to {minLen}.')"
           << "\n";
    python << "    xVal = frames[mask][:minLen]" << "\n";
    python << "    yVal = solverNp[:minLen]" << "\n";
    python << "\n";
    python << "plt.plot(xVal, yVal)" << "\n";
    python << "plt.title('Time used by solver per frame', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.xlabel('Frame', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.ylabel('Solver (ms)', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.grid(True, linestyle=':', linewidth=0.5)" << "\n";
    python << "plt.tight_layout()" << "\n";
    python << "plt.savefig('solver_time_per_frame_py.eps', format='eps', bbox_inches='tight')" << "\n";
    python << "plt.savefig('solver_time_per_frame_py.png', dpi=300, bbox_inches='tight')" << "\n";
    python << "\n";
    python << "\n";
    python << "\n";
    python << "plt.figure(600)" << "\n";
    python << "for i, ylist in enumerate(convergence, start=1):" << "\n";
    python << "    y = np.asarray(ylist, dtype=float)" << "\n";
    python << "\n";
    python << "\n";
    python << "    bad = ~np.isfinite(y) | (y <= 0)" << "\n";
    python << "    if np.any(bad):" << "\n";
    python << "        yClean = y.copy()" << "\n";
    python << "        yClean[bad] = np.nan" << "\n";
    python << "    else:" << "\n";
    python << "        yClean = y" << "\n";
    python << "\n";
    python << "    if yClean.size == 0:" << "\n";
    python << "        continue" << "\n";
    python << "\n";
    python << "    x = np.arange(0, yClean.size + 0)" << "\n";
    python << "    plt.semilogy(x, yClean, alpha=0.8)" << "\n";
    python << "\n";
    python << "plt.title('Convergence', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.xlabel('Solver iteration', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.ylabel('Natural merit function', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.grid(True, which='both', linestyle=':', linewidth=0.5)" << "\n";
    python << "\n";
    python << "plt.savefig('convergence_py.eps', format='eps', bbox_inches='tight')" << "\n";
    python << "plt.savefig('convergence_py.png', dpi=300, bbox_inches='tight')" << "\n";
    python << "\n";
    python << "plt.show()" << "\n";
    python << "\n";
    python << "\n";
    python << "plt.figure(700)" << "\n";
    python << "for i, arr in enumerate(rfactor, start=1):" << "\n";
    python << "    y = toNp(arr)" << "\n";
    python << "    if y.size == 0:" << "\n";
    python << "        continue" << "\n";
    python << "    x = np.arange(0, y.size)" << "\n";
    python << "    plt.plot(x, y, alpha=0.8)" << "\n";
    python << "\n";
    python << "plt.title('R-factor development', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.xlabel('Solver iteration', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.grid(True, linestyle=':', linewidth=0.5)" << "\n";
    python << "plt.tight_layout()" << "\n";
    python << "plt.savefig('rfactors_py.png', dpi=300, bbox_inches='tight')" << "\n";
    python << "plt.show()" << "\n";
    python << "\n";
    python << "maxLen = max((len(c) for c in convergence), default=0)" << "\n";
    python << "if maxLen == 0:" << "\n";
    python << "    raise ValueError('all convergence entries are empty')" << "\n";
    python << "\n";
    python << "data = np.full((len(convergence), maxLen), np.nan, dtype=float)" << "\n";
    python << "for i, ylist in enumerate(convergence):" << "\n";
    python << "    y = np.asarray(ylist, dtype=float)" << "\n";
    python << "    if y.size == 0:" << "\n";
    python << "        continue" << "\n";
    python << "    data[i, :y.size] = y" << "\n";
    python << "    bad = ~np.isfinite(data[i, :y.size]) | (data[i, :y.size] <= 0)" << "\n";
    python << "    if np.any(bad):" << "\n";
    python << "        data[i, :y.size][bad] = np.nan" << "\n";
    python << "\n";
    python << "\n";
    python << "avgCurve = np.nanmean(data, axis=0)" << "\n";
    python << "\n";
    python << "x = np.arange(0, avgCurve.size)" << "\n";
    python << "plt.figure(800)" << "\n";
    python << "plt.semilogy(x, avgCurve, linewidth=0.75)" << "\n";
    python << "plt.title('Average Convergence', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.xlabel('Solver iteration', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.ylabel('Natural merit function', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.grid(True, which='both', linestyle=':', linewidth=0.5)" << "\n";
    python << "\n";
    python << "plt.savefig('convergence_avg.eps', format='eps', bbox_inches='tight')" << "\n";
    python << "plt.savefig('convergence_avg.png', dpi=300, bbox_inches='tight')" << "\n";
    python << "plt.show()" << "\n";
    python << "\n";
    python << "\n";
    python << "convergenceData = np.array(convergence)" << "\n";
    python << "iterationsData = convergenceData.T" << "\n";
    python << "\n";
    python << "q1 = np.percentile(iterationsData, 25, axis=1)" << "\n";
    python << "median = np.percentile(iterationsData, 50, axis=1)" << "\n";
    python << "q3 = np.percentile(iterationsData, 75, axis=1)" << "\n";
    python << "\n";
    python << "plt.figure(figsize=(10, 6))" << "\n";
    python << "plt.plot(q1, label='Q1 (25th Percentile)', color='blue', linestyle='--')" << "\n";
    python << "plt.plot(median, label='Median (50th Percentile)', color='green', linestyle='-')" << "\n";
    python << "plt.plot(q3, label='Q3 (75th Percentile)', color='red', linestyle='--')" << "\n";
    python << "\n";
    python << "plt.title('Quartile Plot for Convergence Across Multiple Frames', fontsize=f_size, fontname=f_type)"
           << "\n";
    python << "plt.xlabel('Solver iteration', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.ylabel('Natural merit function', fontsize=f_size, fontname=f_type)" << "\n";
    python << "plt.yscale('log')" << "\n";
    python << "plt.legend()" << "\n";
    python << "\n";
    python << "plt.savefig('convergence_quartile.eps', format='eps', bbox_inches='tight')" << "\n";
    python << "plt.savefig('convergence_quartile.png', dpi=300, bbox_inches='tight')" << "\n";
    python << "plt.show()" << "\n";
    python.flush();
    python.close();

    logging << "ProxEngine::write_profiling_python(): Done writing profile data..." << newline;

    return true;
}

bool ProxEngine::writeRigidBodiesData(std::string const& filename, unsigned int const& frameNumber)
{
//        typedef typename ProxData::T  T;

    util::Log logging;

    std::string const newline = util::Log::newline();

    std::ofstream python;

        // create directories if necessary
    boost::filesystem::path contact_data_file(filename);

    if (contact_data_file.has_parent_path() && !boost::filesystem::exists(contact_data_file.parent_path())
        && !boost::filesystem::create_directories(contact_data_file.parent_path()))
    {
        util::Log logging;

        logging << "ProxEngine::write_contacts(): Could not create directories for file '" << filename << "'"
                << util::Log::newline();
    }

    python.open(filename.c_str(), std::ios::out);

    auto& bodies = m_data->m_bodies;

    std::vector<std::string> objectNames;
    std::vector<EigenMatrix3<float>> inertiabfs;
    std::vector<float> masses;
    std::vector<EigenQuaternion<float>> orientations;
    std::vector<EigenVector3<float>> positions;
    std::vector<EigenVector3<float>> spin;
    std::vector<EigenVector3<float>> velocities;
    std::vector<std::string> materialNames;
    for (uint32_t i = 0; i < bodies.size(); ++i)
    {
        auto& body = bodies[i];
        inertiabfs.push_back(body.get_inertia_bf());
        masses.push_back(body.get_mass());
        objectNames.push_back("\"" + body.get_name() + "\"");
        orientations.push_back(body.get_orientation());
        positions.push_back(body.get_position());
        spin.push_back(body.get_spin());
        velocities.push_back(body.get_velocity());
        materialNames.push_back("\"" + get_material_name(bodies[i].get_material_idx()) + "\"");
    }

    python << "rigidNames_" << frameNumber << " = " << util::python_write_vector(objectNames) << ";" << "\n";
    python << "inertiabfs_" << frameNumber << " = " << util::python_write_matrix_vector<3>(inertiabfs) << ";" << "\n";
    python << "masses_" << frameNumber << " = " << util::python_write_vector(masses) << ";" << "\n";
    python << "orientations_" << frameNumber << " = " << util::python_write_quaternion(orientations) << ";" << "\n";
    python << "positions_" << frameNumber << " = " << util::python_write_vector(positions) << ";" << "\n";
    python << "spin_" << frameNumber << " = " << util::python_write_vector(spin) << ";" << "\n";
    python << "velocities_" << frameNumber << " = " << util::python_write_vector(velocities) << ";" << "\n";
    python << "materials_" << frameNumber << " = " << util::python_write_vector(materialNames) << ";" << "\n";

        //Same as below, however write to our .py file
    typedef typename ProxData::T T;

    std::vector<T> CX;
    std::vector<T> CY;
    std::vector<T> CZ;
    std::vector<T> NX;
    std::vector<T> NY;
    std::vector<T> NZ;
    std::vector<T> D;
    std::vector<T> AX;
    std::vector<T> AY;
    std::vector<T> AZ;
    std::vector<T> BX;
    std::vector<T> BY;
    std::vector<T> BZ;
    std::vector<unsigned int> A;
    std::vector<unsigned int> B;

    for (const auto& contact : m_data->m_contacts)
    {
        auto [cx, cy, cz] = contact.position;
        auto [nx, ny, nz] = contact.normal;
        auto d = contact.depth;
        auto [ax, ay, az] = (contact.bodyI->get_position());
        auto [bx, by, bz] = (contact.bodyJ->get_position());
        CX.push_back(cx);
        CY.push_back(cy);
        CZ.push_back(cz);
        NX.push_back(nx);
        NY.push_back(ny);
        NZ.push_back(nz);
        D.push_back(d);
        AX.push_back(ax);
        AY.push_back(ay);
        AZ.push_back(az);
        BX.push_back(bx);
        BY.push_back(by);
        BZ.push_back(bz);
        A.push_back(contact.bodyI->get_idx());
        B.push_back(contact.bodyJ->get_idx());
    }

    python << "CX_" << frameNumber << " = " << util::python_write_vector(CX) << ";" << std::endl;
    python << "CY_" << frameNumber << " = " << util::python_write_vector(CY) << ";" << std::endl;
    python << "CZ_" << frameNumber << " = " << util::python_write_vector(CZ) << ";" << std::endl;
    python << "NX_" << frameNumber << " = " << util::python_write_vector(NX) << ";" << std::endl;
    python << "NY_" << frameNumber << " = " << util::python_write_vector(NY) << ";" << std::endl;
    python << "NZ_" << frameNumber << " = " << util::python_write_vector(NZ) << ";" << std::endl;
    python << "D_" << frameNumber << " = " << util::python_write_vector(D) << ";" << std::endl;
    python << "AX_" << frameNumber << " = " << util::python_write_vector(AX) << ";" << std::endl;
    python << "AY_" << frameNumber << " = " << util::python_write_vector(AY) << ";" << std::endl;
    python << "AZ_" << frameNumber << " = " << util::python_write_vector(AZ) << ";" << std::endl;
    python << "BX_" << frameNumber << " = " << util::python_write_vector(BX) << ";" << std::endl;
    python << "BY_" << frameNumber << " = " << util::python_write_vector(BY) << ";" << std::endl;
    python << "BZ_" << frameNumber << " = " << util::python_write_vector(BZ) << ";" << std::endl;
    python << "A_" << frameNumber << " = " << util::python_write_vector(A) << ";" << std::endl;
    python << "B_" << frameNumber << " = " << util::python_write_vector(B) << ";" << std::endl;

    python.flush();
    python.close();

    logging << "ProxEngine::write_contacts(): Done writing rigidbody data..." << newline;

    return true;
}

bool ProxEngine::write_contact_data(std::string const& filename, unsigned int const& frame_number)
{
    typedef typename ProxData::T T;

    util::Log logging;

    std::string const newline = util::Log::newline();

    std::ofstream matlab;

    // create directories if necessary
    boost::filesystem::path contact_data_file(filename);

    if (contact_data_file.has_parent_path() && !boost::filesystem::exists(contact_data_file.parent_path())
        && !boost::filesystem::create_directories(contact_data_file.parent_path()))
    {
        util::Log logging;

        logging << "ProxEngine::write_contacts(): Could not create directories for file '" << filename << "'"
                << util::Log::newline();
    }

    matlab.open(filename.c_str(), std::ios::out);

    if (!matlab.is_open())
    {
        logging << "ProxEngine::write_contacts(): error could not open file = " << filename.c_str()
                << util::Log::newline();

        return false;
    }

    std::vector<T> CX;
    std::vector<T> CY;
    std::vector<T> CZ;
    std::vector<T> NX;
    std::vector<T> NY;
    std::vector<T> NZ;
    std::vector<T> D;
    std::vector<T> AX;
    std::vector<T> AY;
    std::vector<T> AZ;
    std::vector<T> BX;
    std::vector<T> BY;
    std::vector<T> BZ;
    std::vector<unsigned int> A;
    std::vector<unsigned int> B;

    for (const auto& contact : m_data->m_contacts)
    {
        auto [cx, cy, cz] = contact.position;
        auto [nx, ny, nz] = contact.normal;
        auto d = contact.depth;
        auto [ax, ay, az] = (contact.bodyI->get_position());
        auto [bx, by, bz] = (contact.bodyJ->get_position());
        CX.push_back(cx);
        CY.push_back(cy);
        CZ.push_back(cz);
        NX.push_back(nx);
        NY.push_back(ny);
        NZ.push_back(nz);
        D.push_back(d);
        AX.push_back(ax);
        AY.push_back(ay);
        AZ.push_back(az);
        BX.push_back(bx);
        BY.push_back(by);
        BZ.push_back(bz);
        A.push_back(contact.bodyI->get_idx());
        B.push_back(contact.bodyJ->get_idx());
    }

    matlab << "CX_" << frame_number << " = " << util::matlab_write_vector(CX) << ";" << std::endl;
    matlab << "CY_" << frame_number << " = " << util::matlab_write_vector(CY) << ";" << std::endl;
    matlab << "CZ_" << frame_number << " = " << util::matlab_write_vector(CZ) << ";" << std::endl;
    matlab << "NX_" << frame_number << " = " << util::matlab_write_vector(NX) << ";" << std::endl;
    matlab << "NY_" << frame_number << " = " << util::matlab_write_vector(NY) << ";" << std::endl;
    matlab << "NZ_" << frame_number << " = " << util::matlab_write_vector(NZ) << ";" << std::endl;
    matlab << "D_" << frame_number << " = " << util::matlab_write_vector(D) << ";" << std::endl;
    matlab << "AX_" << frame_number << " = " << util::matlab_write_vector(AX) << ";" << std::endl;
    matlab << "AY_" << frame_number << " = " << util::matlab_write_vector(AY) << ";" << std::endl;
    matlab << "AZ_" << frame_number << " = " << util::matlab_write_vector(AZ) << ";" << std::endl;
    matlab << "BX_" << frame_number << " = " << util::matlab_write_vector(BX) << ";" << std::endl;
    matlab << "BY_" << frame_number << " = " << util::matlab_write_vector(BY) << ";" << std::endl;
    matlab << "BZ_" << frame_number << " = " << util::matlab_write_vector(BZ) << ";" << std::endl;
    matlab << "A_" << frame_number << " = " << util::matlab_write_vector(A) << ";" << std::endl;
    matlab << "B_" << frame_number << " = " << util::matlab_write_vector(B) << ";" << std::endl;

    matlab.flush();
    matlab.close();

    logging << "ProxEngine::write_contacts(): Done writing contact data..." << newline;

    return true;
}

} // namespace simulators
