#include "super4pcs/io/io.h"
#include "super4pcs/io/io_ply.h"

#include <Eigen/Geometry>

#ifdef USE_OPENCV
    #include <opencv2/core/core.hpp>
    #include <opencv2/highgui/highgui.hpp>
#endif

#define LINE_BUF_SIZE 100

using namespace std;
using namespace GlobalRegistration;

////////////////////////////////////////////////////////////////////////////////
/// Read
////////////////////////////////////////////////////////////////////////////////
bool
IOManager::ReadObject(const char *name,
           vector<Point3D> &v,
           vector<Eigen::Matrix2f> &tex_coords,
           vector<typename Point3D::VectorType> &normals,
           vector<tripple> &tris,
           vector<std::string> &mtls){
  string filename (name);

  if (filename.length() < 4) return false;

  string ext = filename.substr(filename.size()-3);

  if ( ext.compare ("ply") == 0 )
    return ReadPly (name, v, normals);
  if ( ext.compare ("obj") == 0 )
    return ReadObj (name, v, tex_coords, normals, tris, mtls);
  if ( ext.compare ("ptx") == 0 )
    return ReadPtx (name, v);

  std::cerr << "Unsupported file format" << std::endl;
  return false;
}

bool
IOManager::ReadPly(const char *filename,
                   vector<Point3D> &v,
                   vector<typename Point3D::VectorType> &normals){

  vector<tripple> face;
  unsigned int numOfVertexProperties, numOfVertices, numOfFaces;
  PLYFormat format;
  bool haveColor;
  unsigned int headerSize = readHeader (filename,
                                        numOfVertices,
                                        numOfFaces,
                                        format,
                                        numOfVertexProperties,
                                        haveColor);
  if (haveColor)
    cout << "haveColor" << endl;
  if (headerSize != 0){
    if (format == BINARY_BIG_ENDIAN_1)
      return readBinary1Body (filename, headerSize, numOfVertices,
                       numOfFaces,
                       numOfVertexProperties, haveColor, true, v, normals, face );
    else if (format == BINARY_LITTLE_ENDIAN_1)
      return readBinary1Body (filename, headerSize, numOfVertices,
                       numOfFaces,
                       numOfVertexProperties, haveColor, false, v, normals, face );
    else if (format == ASCII_1)
      return readASCII1Body ( filename, headerSize, numOfVertices,
                       numOfFaces,
                       numOfVertexProperties, haveColor, v, normals, face );
    else{
      cerr << "(PLY) no support for this PLY format" << endl;
      return false;
    }
  }

  return false;
}


bool IOManager::ReadPtx(const char *filename, vector<Point3D> &vertex)
{
    fstream f(filename, ios::in);
    if (!f || f.fail()) {
        cerr << "(PTX) error opening file" << endl;
        return false;
    }


    int numOfVertices;
    int rows, cols;
    char line[LINE_BUF_SIZE];

    {
        f.getline(line,LINE_BUF_SIZE);
        std::stringstream ss(line); ss >> cols;
    }
    {
        f.getline(line,LINE_BUF_SIZE);
        std::stringstream ss(line); ss >> rows;
    }

    numOfVertices = cols*rows;

    // skip matrices declaration
    for(int i=0; i<8; i++) f.getline(line,LINE_BUF_SIZE);

    Point3D ptx;
    float intensity;
    typename Point3D::VectorType rgb;

    vertex.clear();
    vertex.reserve(numOfVertices);


    for (int i = 0; i < numOfVertices && ! f.eof(); i++) {
        f.getline(line,LINE_BUF_SIZE);
        std::stringstream ss(line);

        ss >> ptx.x();
        ss >> ptx.y();
        ss >> ptx.z();
        ss >> intensity;
        ss >> rgb(0);
        ss >> rgb(1);
        ss >> rgb(2);

        ptx.set_rgb(rgb);

        vertex.push_back( ptx );
    }

    f.close();

    return vertex.size() == numOfVertices;
}

bool
IOManager::ReadObj(const char *filename,
                   vector<Point3D> &v,
                   vector<Eigen::Matrix2f> &tex_coords,
                   vector<typename Point3D::VectorType> &normals,
                   vector<tripple> &tris,
                   vector<std::string> &mtls) {
  fstream f(filename, ios::in);
  if (!f || f.fail()) return false;
  char str[1024];
  float x, y, z;
  v.clear();
  tris.clear();
  while (!f.eof()) {
    f.getline(str, 1023);
    char ch[128];
    sscanf(str, "%s %*s", ch);
    if (strcmp(ch, "v") == 0) {
      sscanf(str, "%s %f %f %f", ch, &x, &y, &z);
      v.emplace_back(x, y, z);
      v[v.size() - 1].set_rgb(Point3D::VectorType::Zero());
    } else if (strcmp(ch, "vt") == 0) {
      Eigen::Matrix2f tex_coord;
      sscanf(str, "%s %f %f", ch, &tex_coord.coeffRef(0), &tex_coord.coeffRef(1));
      tex_coords.push_back(tex_coord);
    } else if (strcmp(ch, "vn") == 0) {
      typename Point3D::VectorType normal;
      sscanf(str, "%s %f %f %f", ch, &x, &y, &z);
      normal << x, y, z;
      normals.push_back(normal);
    } else if (strcmp(ch, "f") == 0) {
      tripple triangle;
      if (normals.size() && !tex_coords.size()) {
        sscanf(str, "%s %d//%d %d//%d %d//%d", ch, &(triangle.a),
               &(triangle.n1), &(triangle.b), &(triangle.n2), &(triangle.c),
               &(triangle.n3));
      } else if (normals.size() && tex_coords.size()) {
        sscanf(str, "%s %d/%d/%d %d/%d/%d %d/%d/%d", ch, &(triangle.a),
               &(triangle.t1), &(triangle.n1), &(triangle.b), &(triangle.t2),
               &(triangle.n2), &(triangle.c), &(triangle.t3), &(triangle.n3));
      } else if (!normals.size() && tex_coords.size()) {
        sscanf(str, "%s %d/%d %d/%d %d/%d", ch, &(triangle.a), &(triangle.t1),
               &(triangle.b), &(triangle.t2), &(triangle.c), &(triangle.t3));
      } else if (!normals.size() && !tex_coords.size()) {
        sscanf(str, "%s %d %d %d", ch, &(triangle.a), &(triangle.b),
               &(triangle.c));
      }
      tris.push_back(triangle);
      if (normals.size()) {
        v[triangle.a - 1].set_normal(normals[triangle.n1 - 1]);
        v[triangle.b - 1].set_normal(normals[triangle.n2 - 1]);
        v[triangle.c - 1].set_normal(normals[triangle.n3 - 1]);
      }
    } else if (strcmp(ch, "mtllib") == 0) {
      mtls.push_back(str + 7);
    }
  }
  f.close();

  if(tris.size() == 0){
    // In case we have vertex and normal lists but no face, assign normal to v
    if(v.size() == normals.size()){
      for (size_t i = 0; i < v.size(); ++i)
        v[i].set_normal(normals[i]);
    }
  }else {
    if (! normals.empty()){
      // If we have normals from faces, we must rebuild the normal array to duplicate
      // original normals and get a 1 to 1 correspondances with vertices
      // We assume that the normals have already been sent to vertices
      normals.clear();
      normals.reserve(v.size());

      for (unsigned int i = 0; i!= v.size(); i++)
        normals.push_back(v[i].normal());
    }
  }


  if (mtls.size()) {
    f.open(mtls[0].c_str(), ios::in);
    while (f && !f.fail()) {
      std::string img_name, dummy;
      f >> dummy;
      if (strcmp(dummy.c_str(), "map_Kd") == 0) {
        f >> img_name;
#ifdef USE_OPENCV
        cv::Mat tex = cv::imread(img_name);
        if (!tex.empty()) {
          for (int i = 0; i < tris.size(); ++i) {
            const tripple &t = tris[i];
            Eigen::Matrix2f tc1 = tex_coords[t.t1 - 1];
            Eigen::Matrix2f tc2 = tex_coords[t.t2 - 1];
            Eigen::Matrix2f tc3 = tex_coords[t.t3 - 1];
            if ((tc1.array() < 1.0 && tc1.array() > 1.0 ).all() &&
                (tc2.array() < 1.0 && tc2.array() > 1.0 ).all() &&
                (tc3.array() < 1.0 && tc3.array() > 1.0 ).all()) {

              v[t.a - 1].set_rgb(typename Point3D::VectorType(
                  tex.at<cv::Vec3b>(tc1.coeffRef(1) * tex.rows, tc1.coeffRef(0) * tex.cols)[0],
                  tex.at<cv::Vec3b>(tc1.coeffRef(1) * tex.rows, tc1.coeffRef(0) * tex.cols)[1],
                  tex.at<cv::Vec3b>(tc1.coeffRef(1) * tex.rows, tc1.coeffRef(0) * tex.cols)[2]));

              v[t.b - 1].set_rgb(typename Point3D::VectorType(
                  tex.at<cv::Vec3b>(tc2.coeffRef(1) * tex.rows, tc2.coeffRef(0) * tex.cols)[0],
                  tex.at<cv::Vec3b>(tc2.coeffRef(1) * tex.rows, tc2.coeffRef(0) * tex.cols)[1],
                  tex.at<cv::Vec3b>(tc2.coeffRef(1) * tex.rows, tc2.coeffRef(0) * tex.cols)[2]));

              v[t.c - 1].set_rgb(typename Point3D::VectorType(
                  tex.at<cv::Vec3b>(tc3.coeffRef(1) * tex.rows, tc3.coeffRef(0) * tex.cols)[0],
                  tex.at<cv::Vec3b>(tc3.coeffRef(1) * tex.rows, tc3.coeffRef(0) * tex.cols)[1],
                  tex.at<cv::Vec3b>(tc3.coeffRef(1) * tex.rows, tc3.coeffRef(0) * tex.cols)[2]));
            }
          }
        }
#else
        std::cerr << "OpenCV is required to load material textures. Skipping "
                  << img_name.c_str()
                  << std::endl;
#endif
      }
    }
  }
  f.close();

  if (v.size() == 0) return false;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Write
////////////////////////////////////////////////////////////////////////////////

bool
IOManager::WriteObject(const char *name,
                       const vector<Point3D> &v,
                       const vector<Eigen::Matrix2f> &tex_coords,
                       const vector<typename Point3D::VectorType> &normals,
                       const vector<tripple> &tris,
                       const vector<std::string> &mtls) {
  string filename (name);
  string ext = filename.substr(filename.size()-3);

  bool haveExt = filename.at(filename.size()-4) == '.';

  if (tris.size() == 0){
    return WritePly(haveExt ?
                    filename.substr(0,filename.size()-3).append("ply") :
                    filename.append(".ply"),
                    v, normals);
  }
  else{
    return WriteObj(haveExt ?
                    filename.substr(0,filename.size()-3).append("obj") :
                    filename.append(".obj"),
                    v,
                    tex_coords,
                    normals,
                    tris,
                    mtls);
  }

}

bool IOManager::WriteMatrix(
        const string &name,
        const Eigen::Ref<const Eigen::Matrix<double, 4, 4> >& mat,
        IOManager::MATRIX_MODE mode)
{
    std::ofstream sstr;
    sstr.open(name, std::ofstream::out | std::ofstream::trunc);

    bool status = false;

    switch (mode) {
    case POLYWORKS:
        formatPolyworksMatrix(mat, sstr);
        status = true;
        break;
    default:
        break;
    }

    sstr.close();

    return status;
}

bool
IOManager::WritePly(string filename,
                    const vector<Point3D> &v,
                    const vector<typename Point3D::VectorType> &normals)
{
  std::ofstream plyFile;
  plyFile.open (filename.c_str(), std::ios::out |  std::ios::trunc | std::ios::binary);
  if (! plyFile.is_open()){
    std::cerr << "Cannot open file to write!" << std::endl;
    return false;
  }

  bool useNormals = normals.size() == v.size();
  // we check if we have colors by looking if the first rgb vector is void
  bool useColors = false;
  for (unsigned int i = 0; i!=v.size(); i++){
    if (v[i].hasColor()){
      useColors = true;
      break;
    }
  }

  plyFile.imbue(std::locale::classic());

  // Write Header
  plyFile << "ply" << std::endl;
  plyFile << "format binary_little_endian 1.0" << std::endl;
  plyFile << "comment Super4PCS output file" << std::endl;
  plyFile << "element vertex " << v.size() << std::endl;
  plyFile << "property float x" << std::endl;
  plyFile << "property float y" << std::endl;
  plyFile << "property float z" << std::endl;

  if(useNormals) {
    plyFile << "property float nx" << std::endl;
    plyFile << "property float ny" << std::endl;
    plyFile << "property float nz" << std::endl;
  }

  if(useColors) {
    plyFile << "property uchar red" << std::endl;
    plyFile << "property uchar green" << std::endl;
    plyFile << "property uchar blue" << std::endl;
  }

  plyFile << "end_header" << std::endl;

  // Read all elements in data, correct their depth and print them in the file
  char tmpChar;
  float tmpFloat;
  for (unsigned int i = 0; i!=v.size(); i++){
    tmpFloat = v[i].x();
    plyFile.write(reinterpret_cast<const char*>(&tmpFloat),sizeof(float));
    tmpFloat = v[i].y();
    plyFile.write(reinterpret_cast<const char*>(&tmpFloat),sizeof(float));
    tmpFloat = v[i].z();
    plyFile.write(reinterpret_cast<const char*>(&tmpFloat),sizeof(float));

    if (useNormals){
      plyFile.write(reinterpret_cast<const char*>(&normals[i](0)),sizeof(float));
      plyFile.write(reinterpret_cast<const char*>(&normals[i](1)),sizeof(float));
      plyFile.write(reinterpret_cast<const char*>(&normals[i](2)),sizeof(float));
    }

    if (useColors){
      tmpChar = v[i].rgb()[0];
      plyFile.write(reinterpret_cast<const char*>(&tmpChar),sizeof(char));
      tmpChar = v[i].rgb()[1];
      plyFile.write(reinterpret_cast<const char*>(&tmpChar),sizeof(char));
      tmpChar = v[i].rgb()[2];
      plyFile.write(reinterpret_cast<const char*>(&tmpChar),sizeof(char));
    }
  }

  plyFile.close();

  return true;
}

bool
IOManager::WriteObj(string filename, const vector<Point3D> &v,
                    const vector<Eigen::Matrix2f> &tex_coords,
                    const vector<typename Point3D::VectorType> &normals,
                    const vector<tripple> &tris,
                    const vector<std::string> &mtls) {
  fstream f(filename.c_str(), ios::out);
  if (!f || f.fail()) return false;
  size_t i;

  //normals.clear();

  for (i = 0; i < mtls.size(); ++i) {
    f << "mtllib " << mtls[i] << endl;
  }

  for (i = 0; i < v.size(); ++i) {
    f << "v "
      << v[i].x() << " " << v[i].y() << " " << v[i].z() << " ";

    if (v[i].rgb()[0] != 0)
      f << v[i].rgb()[0] << " " << v[i].rgb()[1] << " " << v[i].rgb()[2];

    f << endl;
  }

  for (i = 0; i < normals.size(); ++i) {
    f << "vn " << normals[i](0) << " " << normals[i](1) << " " << normals[i](2)
      << endl;
  }

  for (i = 0; i < tex_coords.size(); ++i) {
    f << "vt " << tex_coords[i].coeffRef(0) << " " << tex_coords[i].coeffRef(1) << endl;
  }

  for (i = 0; i < tris.size(); ++i) {
    if (!normals.size() && !tex_coords.size())
      f << "f " << tris[i].a << " " << tris[i].b << " " << tris[i].c << endl;
    else if (tex_coords.size())
      f << "f " << tris[i].a << "/" << tris[i].t1 << " " << tris[i].b << "/"
        << tris[i].t2 << " " << tris[i].c << "/" << tris[i].t3 << endl;
    else
        f << "f " << tris[i].a << "/" << tris[i].n1 << " " << tris[i].b << "/"
        << tris[i].n2 << " " << tris[i].c << "/" << tris[i].n3 << endl;
  }

  f.close();

  return true;
}



std::ofstream &
IOManager::formatPolyworksMatrix(
        const Eigen::Ref<const Eigen::Matrix<double, 4, 4> > &mat,
        std::ofstream &sstr) {

    auto formatValue = [](const double& v){
        return (v >= 0.
                ? std::string(" ") + std::to_string(v)
                : std::to_string(v));
    };

    sstr << "VERSION\t=\t1\n";
    sstr << "MATRIX\t=\n";
    for (int j = 0; j!= 4; ++j){
        sstr << formatValue(mat(j, 0)) << "  "
             << formatValue(mat(j, 1)) << "  "
             << formatValue(mat(j, 2)) << "  "
             << formatValue(mat(j, 3)) << "\n";
    }

    return sstr;
}




