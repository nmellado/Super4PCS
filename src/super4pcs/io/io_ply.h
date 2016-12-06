#ifndef IO_PLY_H
#define IO_PLY_H

typedef enum {
    BINARY_BIG_ENDIAN_1,
    BINARY_LITTLE_ENDIAN_1,
    ASCII_1
} PLYFormat;

static const unsigned int MAX_COMMENT_SIZE = 256;


unsigned int
readHeader ( const char *filename,
             unsigned int & numOfVertices,
             unsigned int & numOfFaces,
             PLYFormat & format,
             unsigned int & numOfVertexProperties,
             bool         & haveColor )
{
    // use binary mode to preserve end line symbols on linux and windows
    ifstream in (filename, std::ios_base::in | std::ios_base::binary);
    if (!in){
        cerr << "(PLY) error opening file" << endl;
        return 0;
    }

    numOfVertexProperties = 0;
    numOfVertices = 0;
    numOfFaces = 0;
    haveColor = false;

    string current, currentelement;
    in >> current;
    if (current != "ply"){
        cerr << "(PLY) not a PLY file" << endl;
        return 0;
    }
    in >> current;
    int lid = 0;
    while (current != "end_header") {
        if (current == "format") {
            in >> current;
            if (current == "binary_big_endian") {
                in >> current;
                if (current == "1.0")
                    format = BINARY_BIG_ENDIAN_1;
                else{
                    cerr << "(PLY) error parsing header - bad binary big endian version" << endl;
                    return 0;
                }
            } else if (current == "binary_little_endian") {
                in >> current;
                if (current == "1.0")
                    format = BINARY_LITTLE_ENDIAN_1;
                else{
                    cerr << "(PLY) error parsing header - bad binary little endian version" << endl;
                    return 0;
                }
            } else if (current == "ascii") {
                in >> current;
                if (current == "1.0")
                    format = ASCII_1;
                else{
                    cerr << "(PLY) error parsing header - bad ascii version" << endl;
                    return 0;
                }
            } else {
                cerr << "(PLY) error parsing header (format)" << endl;
                return 0;
            }
        } else if (current == "element") {
            in >> current;
            if (current == "vertex"){
                currentelement = current;
                in >> numOfVertices;
            }
            else if (current == "face"){
                currentelement = current;
                in >> numOfFaces;
            }
            else{
                cerr << "(PLY) ignoring unknown element " << current << endl;
                currentelement = "";
            }
        } else if (currentelement != "" && current == "property") {
            in >> current;
            if (current == "float") {
                numOfVertexProperties++;
                in >> current;
            }
            else if (current == "uchar") { // color
                numOfVertexProperties++;
                haveColor = true;
                in >> current;
            }
            else if (current == "list") {
                in >> current;
                in >> current;
                in >> current;
            } else {
                cerr << "(PLY) error parsing header (property)" << endl;
                return 0;
            }
        } else if ( (current == "comment") || (current.find("obj_info") != std::string::npos) ) {
            char comment[MAX_COMMENT_SIZE];
            in.getline (comment, MAX_COMMENT_SIZE);
        } else {
            ;
            //cerr << "(PLY) ignoring line: " << current << endl;
            //return 0;
        }
        in >> current;
        lid++;
    }

    unsigned int headerSize = in.tellg ();
    in.close ();
    return headerSize+1;
}


template <class T>
void
bigLittleEndianSwap (T * v, unsigned int numOfElements)
{
    char * tmp = (char*)v;
    for (unsigned int j = 0; j < numOfElements; j++){
        unsigned int offset = 4*j;
        char c = tmp[offset];
        tmp[offset] =  tmp[offset+3];
        tmp[offset+3] = c;
        c = tmp[offset+1];
        tmp[offset+1] = tmp[offset+2];
        tmp[offset+2] = c;
    }
}


bool
readBinary1Body (const std::string & filename,
                 unsigned int headerSize,
                 unsigned int numOfVertices,
                 unsigned int numOfFaces,
                 unsigned int numOfVertexProperties,
                 bool         haveColor,
                 bool bigEndian,
                 vector<Point3D>& vertex,
                 vector<cv::Point3f>& normal,
                 vector<tripple>& face )
{
    //size_t count;

    FILE * in = fopen (filename.c_str (), "r");
    if (!in){
        cerr << "(PLY) error opening file" << endl;
        return false;
    }

    for (unsigned int i = 0; i < headerSize; i++) {
        char c;
        /*count = */fread (&c, 1, 1, in);
    }

    // *****************
    // Reading geometry.
    // *****************
    cv::Point3f n;
    cv::Vec3f rgb;
    float * v = new float[numOfVertexProperties];
    uchar rgb_buff [4];

    for (unsigned int i = 0; i < numOfVertices && !feof (in); i++) {
        if (numOfVertexProperties==10){
            fread (v, 4, 6, in);
            fread (rgb_buff, sizeof(uchar), 4, in);
        }else if (numOfVertexProperties==9){
            fread (v, 4, 6, in);
            fread (rgb_buff, sizeof(uchar), 3, in);
        }else if (numOfVertexProperties==6 && haveColor){
            fread (v, 4, 3, in);
            fread (rgb_buff, sizeof(uchar), 3, in);
        }else if (numOfVertexProperties==7 ){
            fread (v, 4, 3, in);
            fread (rgb_buff, sizeof(uchar), 4, in);
        }
        else
            fread (v, 4, numOfVertexProperties, in);
        if (bigEndian == true)
            bigLittleEndianSwap (v, numOfVertexProperties);
        vertex.push_back( Point3D(v[0],v[1],v[2]) );

        if (numOfVertexProperties == 6){
            if (haveColor){
                rgb[0] = rgb_buff[0];
                rgb[1] = rgb_buff[1];
                rgb[2] = rgb_buff[2];
                vertex.back().set_rgb(rgb);
            }else{
                n.x = v[3];
                n.y = v[4];
                n.z = v[5];
                normal.push_back (n);
                vertex.back().set_normal(n);
            }
        }else if (numOfVertexProperties == 7){
            rgb[0] = rgb_buff[0];
            rgb[1] = rgb_buff[1];
            rgb[2] = rgb_buff[2];
            vertex.back().set_rgb(rgb);
        }else if (numOfVertexProperties == 9 || numOfVertexProperties == 10){
            n.x = v[3];
            n.y = v[4];
            n.z = v[5];
            rgb[0] = rgb_buff[0];
            rgb[1] = rgb_buff[1];
            rgb[2] = rgb_buff[2];
            normal.push_back (n);
            vertex.back().set_normal(n);
            vertex.back().set_rgb(rgb);
        }
    }
    delete [] v;

    if (numOfFaces != 0){
        if (feof (in)){
            cerr << "(PLY) incomplete file" << endl;
            return false;
        }

        // *****************
        // Reading topology.
        // *****************
        for (unsigned int i = 0; i < numOfFaces && !feof (in); i++) {
            unsigned int f[4];
            char polygonSize;
            /*count = */fread (&polygonSize, 1, 1, in);
            /*count = */fread (f, 4, 3, in);
            if (bigEndian == true)
                bigLittleEndianSwap (f, 3);
            face.push_back(tripple(f[0],f[1],f[2]));
        }
    }

    return true;
}


bool
readASCII1Body (const std::string & filename,
                unsigned int headerSize,
                unsigned int numOfVertices,
                unsigned int numOfFaces,
                unsigned int numOfVertexProperties,
                bool         haveColor,
                vector<Point3D>& vertex,
                vector<cv::Point3f>& normal,
                vector<tripple>& face )
{

    FILE * in = fopen (filename.c_str (), "r");
    if (!in){
        cerr << "(PLY) error opening file" << endl;
        return false;
    }

    for (unsigned int i = 0; i < headerSize; i++) {
        char c;
        /*count = */fread (&c, 1, 1, in);
    }

    // *****************
    // Reading geometry.
    // *****************
    cv::Point3f n;
    cv::Vec3f rgb;
    unsigned int rgb_buff [4];
    for (unsigned int i = 0; i < numOfVertices && !feof (in); i++) {
        std::vector<float> v(numOfVertexProperties);

        if (numOfVertexProperties==10){
            for (unsigned int j = 0;  j < 6;  j++)
                fscanf (in, "%f", &v[j]);
            for (unsigned int j = 0;  j < 4;  j++)
                fscanf (in, "%i", &rgb_buff[j]);
        }
        else if (numOfVertexProperties==9){
            for (unsigned int j = 0;  j < 6;  j++)
                fscanf (in, "%f", &v[j]);
            for (unsigned int j = 0;  j < 3;  j++)
                fscanf (in, "%i", &rgb_buff[j]);
        }else if (numOfVertexProperties==6 && haveColor){
            for (unsigned int j = 0;  j < 3;  j++)
                fscanf (in, "%f", &v[j]);
            for (unsigned int j = 0;  j < 3;  j++)
                fscanf (in, "%i", &rgb_buff[j]);
        }else if (numOfVertexProperties==7){
            for (unsigned int j = 0;  j < 3;  j++)
                fscanf (in, "%f", &v[j]);
            for (unsigned int j = 0;  j < 4;  j++)
                fscanf (in, "%i", &rgb_buff[j]);
        }
        else
            for (unsigned int j = 0;  j < numOfVertexProperties;  j++)
                fscanf (in, "%f", &v[j]);
        
        vertex.push_back( Point3D(v[0],v[1],v[2]) );

        if (numOfVertexProperties == 6){
            if (haveColor){
                rgb[0] = float(rgb_buff[0]);
                rgb[1] = float(rgb_buff[1]);
                rgb[2] = float(rgb_buff[2]);
                vertex.back().set_rgb(rgb);

            }else{
                n.x = v[3];
                n.y = v[4];
                n.z = v[5];
                normal.push_back (n);
                vertex.back().set_normal(n);
            }
        }else if (numOfVertexProperties == 7){
            rgb[0] = float(rgb_buff[0]);
            rgb[1] = float(rgb_buff[1]);
            rgb[2] = float(rgb_buff[2]);
            vertex.back().set_rgb(rgb);
        }else if (numOfVertexProperties == 9 || numOfVertexProperties == 10){
            n.x = v[3];
            n.y = v[4];
            n.z = v[5];
            rgb[0] = float(rgb_buff[0]);
            rgb[1] = float(rgb_buff[1]);
            rgb[2] = float(rgb_buff[2]);
            normal.push_back (n);
            vertex.back().set_normal(n);
            vertex.back().set_rgb(rgb);
        }
    }

    if (numOfFaces != 0){
        if (feof (in)){
            cerr << "(PLY) incomplete file" << endl;
            return false;
        }

        // *****************
        // Reading topology.
        // *****************
        for (unsigned int i = 0; i < numOfFaces && !feof (in); i++) {
            int f[3];
            int polygonSize;
            /*count = */fscanf (in, "%d %d %d %d", &polygonSize, &f[0], &f[1], &f[2]);
            face.push_back(tripple(f[0],f[1],f[2]));
        }
    }

    return true;
}

#endif //IO_PLY_H
