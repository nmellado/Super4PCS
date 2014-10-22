#include "io.h"
#include "io_ply.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
/// Read
////////////////////////////////////////////////////////////////////////////////
bool 
IOManager::ReadObject(
           const char *name,
           vector<Point3D> &v, 
           vector<cv::Point2f> &tex_coords,
           vector<cv::Point3f> &normals, 
           vector<tripple> &tris,
           vector<std::string> &mtls){
  string filename (name);
  string ext = filename.substr(filename.size()-3);
  
  if ( ext.compare ("ply") == 0 )
    return ReadPly (name, v, normals);
  if ( ext.compare ("obj") == 0 )
    return ReadObj (name, v, tex_coords, normals, tris, mtls);
  
  std::cerr << "Unsupported file format" << std::endl;  
  return false;
}

bool 
IOManager::ReadPly(const char *filename,
                   vector<Point3D> &v, 
                   vector<cv::Point3f> &normals){
  
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

bool 
IOManager::ReadObj(const char *filename,
                   vector<Point3D> &v, 
                   vector<cv::Point2f> &tex_coords,
                   vector<cv::Point3f> &normals, 
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
      v.push_back(Point3D(x, y, z));
      v[v.size() - 1].set_rgb(cv::Vec3f(0, 0, 0));
    } else if (strcmp(ch, "vt") == 0) {
      cv::Point2f tex_coord;
      sscanf(str, "%s %f %f", ch, &tex_coord.x, &tex_coord.y);
      tex_coords.push_back(tex_coord);
    } else if (strcmp(ch, "vn") == 0) {
      cv::Point3f normal;
      sscanf(str, "%s %f %f %f", ch, &(normal.x), &(normal.y), &(normal.z));
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
      for (int i = 0; i < v.size(); ++i) 
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
        cv::Mat tex = cv::imread(img_name);
        if (!tex.empty()) {
          for (int i = 0; i < tris.size(); ++i) {
            const tripple &t = tris[i];
            cv::Point2f tc1 = tex_coords[t.t1 - 1];
            cv::Point2f tc2 = tex_coords[t.t2 - 1];
            cv::Point2f tc3 = tex_coords[t.t3 - 1];
            if (tc1.x < 1.0 && tc1.x > 0 && tc1.y < 1.0 && tc1.y > 0 &&
                tc2.x < 1.0 && tc2.x > 0 && tc2.y < 1.0 && tc2.y > 0 &&
                tc3.x < 1.0 && tc3.x > 0 && tc3.y < 1.0 && tc3.y > 0) {
                
              v[t.a - 1].set_rgb(cv::Vec3f(
                  static_cast<float>(
                      tex.at<cv::Vec3b>(tc1.y * tex.rows, tc1.x * tex.cols)[0]),
                  static_cast<float>(
                      tex.at<cv::Vec3b>(tc1.y * tex.rows, tc1.x * tex.cols)[1]),
                  static_cast<float>(tex.at<cv::Vec3b>(tc1.y * tex.rows,
                                                       tc1.x * tex.cols)[2])));

              v[t.b - 1].set_rgb(cv::Vec3f(
                  static_cast<float>(
                      tex.at<cv::Vec3b>(tc2.y * tex.rows, tc2.x * tex.cols)[0]),
                  static_cast<float>(
                      tex.at<cv::Vec3b>(tc2.y * tex.rows, tc2.x * tex.cols)[1]),
                  static_cast<float>(tex.at<cv::Vec3b>(tc2.y * tex.rows,
                                                       tc2.x * tex.cols)[2])));

              v[t.c - 1].set_rgb(cv::Vec3f(
                  static_cast<float>(
                      tex.at<cv::Vec3b>(tc3.y * tex.rows, tc3.x * tex.cols)[0]),
                  static_cast<float>(
                      tex.at<cv::Vec3b>(tc3.y * tex.rows, tc3.x * tex.cols)[1]),
                  static_cast<float>(tex.at<cv::Vec3b>(tc3.y * tex.rows,
                                                       tc3.x * tex.cols)[2])));
            }
          }
        }
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
                       vector<Point3D> &v,
                       vector<cv::Point2f> &tex_coords, 
                       vector<cv::Point3f> &normals,
                       vector<tripple> &tris, 
                       vector<std::string> &mtls) {
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
bool
IOManager::WritePly(string filename, 
                    vector<Point3D> &v, 
                    vector<cv::Point3f> &normals)
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
    if (cv::norm(v[i].rgb())){
      useColors = true;
      break;
    }
  }
  
  plyFile.imbue(std::locale::classic());
  
  // Write Header
  plyFile << "ply" << std::endl;
  plyFile << "format binary_little_endian 1.0" << std::endl;
  plyFile << "comment Kinect depth data" << std::endl;
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
  char tmpColor;
  for (unsigned int i = 0; i!=v.size(); i++){
    plyFile.write(reinterpret_cast<char*>(&v[i].x),sizeof(float));
    plyFile.write(reinterpret_cast<char*>(&v[i].y),sizeof(float));
    plyFile.write(reinterpret_cast<char*>(&v[i].z),sizeof(float));

    if (useNormals){
      plyFile.write(reinterpret_cast<char*>(&normals[i].x),sizeof(float));
      plyFile.write(reinterpret_cast<char*>(&normals[i].y),sizeof(float));
      plyFile.write(reinterpret_cast<char*>(&normals[i].z),sizeof(float));
    }

    if (useColors){
      tmpColor = v[i].rgb()[0];
      plyFile.write(reinterpret_cast<char*>(&tmpColor),sizeof(char));
      tmpColor = v[i].rgb()[1];
      plyFile.write(reinterpret_cast<char*>(&tmpColor),sizeof(char));
      tmpColor = v[i].rgb()[2];
      plyFile.write(reinterpret_cast<char*>(&tmpColor),sizeof(char));
    }
  }
    
  plyFile.close();
  
  cout << "Merged object was written to " << filename.c_str() << endl;
  return true;
}

bool 
IOManager::WriteObj(string filename, vector<Point3D> &v,
                    vector<cv::Point2f> &tex_coords, 
                    vector<cv::Point3f> &normals,
                    vector<tripple> &tris, 
                    vector<std::string> &mtls) {
  fstream f(filename.c_str(), ios::out);
  if (!f || f.fail()) return false;
  int i;

  normals.clear();

  for (i = 0; i < mtls.size(); ++i) {
    f << "mtllib " << mtls[i] << endl;
  }

  for (i = 0; i < v.size(); ++i) {
    f << "v " 
      << v[i].x << " " << v[i].y << " " << v[i].z << " ";
      
    if (v[i].rgb()[0] != 0)
      f << v[i].rgb()[0] << " " << v[i].rgb()[1] << " " << v[i].rgb()[2];
    
    f << endl;
  }

  for (i = 0; i < normals.size(); ++i) {
    f << "vn " << normals[i].x << " " << normals[i].y << " " << normals[i].z
      << endl;
  }

  for (i = 0; i < tex_coords.size(); ++i) {
    f << "vt " << tex_coords[i].x << " " << tex_coords[i].y << endl;
  }

  for (i = 0; i < tris.size(); ++i) {
    if (!normals.size() && !tex_coords.size())
      f << "f " << tris[i].a << " " << tris[i].b << " " << tris[i].c << endl;
    else if (tex_coords.size())
      f << "f " << tris[i].a << "/" << tris[i].t1 << " " << tris[i].b << "/"
        << tris[i].t2 << " " << tris[i].c << "/" << tris[i].t3 << endl;
  }

  f.close();
  
  cout << "Merged object was written to " << filename.c_str() << endl;
  return true;
}






