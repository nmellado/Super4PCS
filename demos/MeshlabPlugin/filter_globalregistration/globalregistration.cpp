/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#include "globalregistration.h"
#include "super4pcs/algorithms/super4pcs.h"
#include "super4pcs/algorithms/4pcs.h"
#include <QtScript>

// Constructor usually performs only two simple tasks of filling the two lists
//  - typeList: with all the possible id of the filtering actions
//  - actionList with the corresponding actions. If you want to add icons to your filtering actions you can do here by construction the QActions accordingly

GlobalRegistrationPlugin::GlobalRegistrationPlugin()
{
    typeList << FP_GLOBAL_REGISTRATION;

  foreach(FilterIDType tt , types())
      actionList << new QAction(filterName(tt), this);
}

// ST() must return the very short string describing each filtering action
// (this string is used also to define the menu entry)
QString GlobalRegistrationPlugin::filterName(FilterIDType filterId) const
{
  switch(filterId) {
        case FP_GLOBAL_REGISTRATION :  return QString("Global registration");
        default : assert(0);
    }
  return QString();
}

// Info() must return the longer string describing each filtering action
// (this string is used in the About plugin dialog)
 QString GlobalRegistrationPlugin::filterInfo(FilterIDType filterId) const
{
  switch(filterId) {
        case FP_GLOBAL_REGISTRATION :  return QString("Compute the rigid transforation aligning two 3d objets.");
        default : assert(0);
    }
    return QString("Unknown Filter");
}

// The FilterClass describes in which generic class of filters it fits.
// This choice affect the submenu in which each filter will be placed
// More than a single class can be choosen.
GlobalRegistrationPlugin::FilterClass GlobalRegistrationPlugin::getClass(QAction *a)
{
  switch(ID(a))
    {
        case FP_GLOBAL_REGISTRATION :  return MeshFilterInterface::PointSet;
        default : assert(0);
    }
    return MeshFilterInterface::Generic;
}

// This function define the needed parameters for each filter. Return true if the filter has some parameters
// it is called every time, so you can set the default value of parameters according to the mesh
// For each parameter you need to define,
// - the name of the parameter,
// - the string shown in the dialog
// - the default value
// - a possibly long string describing the meaning of that parameter (shown as a popup help in the dialog)
void GlobalRegistrationPlugin::initParameterSet(QAction *action,MeshDocument &md, RichParameterSet & parlst)
{

     switch(ID(action))	 {
        case FP_GLOBAL_REGISTRATION :

         parlst.addParam(new RichMesh ("refMesh",md.mm(),&md, "Reference Mesh",	"Reference point-cloud or mesh"));
         parlst.addParam(new RichMesh ("targetMesh",md.mm(),&md, "Target Mesh",	"Point-cloud or mesh to be aligned to the reference"));
         parlst.addParam(new RichAbsPerc("overlap", 50, 0, 100, "Overlap Ratio", "Overlap ratio between the two clouds (command line option: -o)"));
         parlst.addParam(new RichFloat("delta",   0.1, "Registration tolerance", "Tolerance value for the congruent set exploration and LCP computation (command line option: -d)"));
         parlst.addParam(new RichInt("nbSamples", 200, "Number of samples", "Number of samples used in each mesh (command line option: -n)"));
         parlst.addParam(new RichFloat("norm_diff", -1, "Filter: difference of normal (degrees)", "Allowed difference of normals allowed between corresponding pairs of points(command line option: -a)"));
         parlst.addParam(new RichFloat("color_diff", -1, "Filter: difference color", "Allowed difference of colors allowed between corresponding pairs of points(command line option: -c)"));
         parlst.addParam(new RichInt("max_time_seconds", 10000, "Max. Computation time, in seconds", "Stop the computation before the end of the exploration (command line option: -t)"));
         parlst.addParam(new RichBool("useSuper4PCS", true, "Use Super4PCS", "When disable, use 4PCS algorithm (command line option: -x"));

         break;
     default : assert(0);
    }
}

struct RealTimeTransformVisitor {
    using MatrixType = typename GlobalRegistration::Match4PCSBase::MatrixType;
    CMeshO* mesh = nullptr;
    GlobalRegistrationPlugin* plugin;
    inline void operator() (
            float /*fraction*/,
            float best_LCP,
            Eigen::Ref<MatrixType> mat) {
        plugin->Log("Found new configuration. LCP = %f", best_LCP);

        mesh->Tr.FromEigenMatrix(mat);
    }
    constexpr bool needsGlobalTransformation() { return true; }
};

struct TransformVisitor {
    using MatrixType = typename GlobalRegistration::Match4PCSBase::MatrixType;
    CMeshO* mesh = nullptr;
    GlobalRegistrationPlugin* plugin;
    inline void operator() (
            float /*fraction*/,
            float best_LCP,
            Eigen::Ref<MatrixType> /*mat*/) const {
        plugin->Log("Found new configuration. LCP = %f", best_LCP);
    }
    constexpr bool needsGlobalTransformation() { return false; }
};

// The Real Core Function doing the actual mesh processing.
// Move Vertex of a random quantity
bool GlobalRegistrationPlugin::applyFilter(QAction */*filter*/,
                                           MeshDocument &/*md*/,
                                           RichParameterSet & par,
                                           vcg::CallBackPos */*cb*/)
{

    MeshModel *mmref = par.getMesh("refMesh");
    MeshModel *mmtrg = par.getMesh("targetMesh");
    CMeshO *refMesh=&mmref->cm;
    CMeshO *trgMesh=&mmtrg->cm;

    using Sampler = GlobalRegistration::Sampling::UniformDistSampler;

//    Log("Initializing Super4PCS. Delta=%f, overlap=%f", delta, overlap);

    GlobalRegistration::Match4PCSOptions opt;
    opt.configureOverlap(par.getAbsPerc("overlap")/100.f);
    opt.delta                 = par.getFloat("delta");
    opt.sample_size           = par.getInt("nbSamples");
    opt.max_normal_difference = par.getFloat("norm_diff");
    opt.max_color_distance    = par.getFloat("color_diff");
    opt.max_time_seconds      = par.getInt("max_time_seconds");

    bool useSuper4PCS         = par.getBool("useSuper4PCS");


    GlobalRegistration::Utils::Logger logger (GlobalRegistration::Utils::LogLevel::NoLog);
    GlobalRegistration::Match4PCSBase* matcher = nullptr;
    Sampler sampler;

    if (useSuper4PCS)
        matcher = new GlobalRegistration::MatchSuper4PCS (opt, logger);
    else
        matcher = new GlobalRegistration::Match4PCS (opt, logger);

    GlobalRegistration::Match4PCSBase::MatrixType mat;
    std::vector<GlobalRegistration::Point3D> set1, set2;

    // init Super4PCS point cloud internal structure
    auto fillPointSet = [] (const CMeshO& m, std::vector<GlobalRegistration::Point3D>& out) {
        using GlobalRegistration::Point3D;
        Point3D p;
        out.clear();
        out.reserve(m.vert.size());

        // TODO: copy other point-wise information, if any
        for(size_t i = 0; i< m.vert.size(); i++){
            const auto& vertex = m.vert[i];
            vertex.P().ToEigenVector(p.pos());
            out.push_back(p);
        }
    };
    fillPointSet(*refMesh, set1);
    fillPointSet(*trgMesh, set2);

    // run
    TransformVisitor v;
    v.mesh = trgMesh;
    v.plugin = this;

    float score = matcher->ComputeTransformation(set1, &set2, mat, sampler, v);
    Log("Final LCP = %f", score);
    v.mesh->Tr.FromEigenMatrix(mat);

    delete matcher;

    return true;
}

MESHLAB_PLUGIN_NAME_EXPORTER(GlobalRegistrationPlugin)
