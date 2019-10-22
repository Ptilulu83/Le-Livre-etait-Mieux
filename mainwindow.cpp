#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>

/* **** début de la partie à compléter **** */


float MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    FaceHandle fh = _mesh->face_handle ( faceID );

    std::vector<VertexHandle> vertexes;

    MyMesh::FaceVertexIter fh_v = _mesh->fv_iter(fh);
    for(; fh_v.is_valid(); ++fh_v)
        vertexes.push_back ( *fh_v );

    float valuesAB[3];
    valuesAB[0] = _mesh->point(vertexes[1])[0] - _mesh->point(vertexes[0])[0];
    valuesAB[1] = _mesh->point(vertexes[1])[1] - _mesh->point(vertexes[0])[1];
    valuesAB[2] = _mesh->point(vertexes[1])[2] - _mesh->point(vertexes[0])[2];

    float valuesAC[3];
    valuesAC[0] = _mesh->point(vertexes[2])[0] - _mesh->point(vertexes[0])[0];
    valuesAC[1] = _mesh->point(vertexes[2])[1] - _mesh->point(vertexes[0])[1];
    valuesAC[2] = _mesh->point(vertexes[2])[2] - _mesh->point(vertexes[0])[2];

    VectorT<float, 3> vectorAB(valuesAB);
    VectorT<float, 3> vectorAC(valuesAC);

    VectorT<float, 3> product = vectorAB * vectorAC;
    float norm = product.norm();

    return norm / 2.0f;
}

float MainWindow::baryArea(MyMesh* _mesh, int vertID){
    float baryArea = 0;

    VertexHandle vh = _mesh->vertex_handle ( vertID );
    MyMesh::VertexFaceIter vf = _mesh->vf_iter ( vh );
    for ( ; vf.is_valid ( ) ; ++vf ) {
        FaceHandle current = *vf;
        baryArea += faceArea ( _mesh , current.idx( ) );
    }
    return baryArea / 3.0f;
}

float MainWindow::angleFF(MyMesh* _mesh, int faceID0,  int faceID1, int vertID0, int vertID1)
{
    int sign = 0;

    VertexHandle vh0 = _mesh->vertex_handle ( vertID0 );
    FaceHandle fh0 = _mesh->face_handle ( faceID0 );

    MyMesh::FaceVertexCWIter fh_cwv = _mesh->fv_cwiter ( fh0 );
    while ( fh_cwv.is_valid ( ) && *fh_cwv != vh0 ) ++fh_cwv;

    VertexHandle next = *++fh_cwv;

    if ( next.idx ( ) == vertID1 ) sign = -1;
    else sign = 1;

    OpenMesh::Vec3f normal0 (_mesh->normal ( fh0 ) );
    OpenMesh::Vec3f normal1 (_mesh->normal ( _mesh->face_handle ( faceID1 ) ) );

    float scalar = normal0 | normal1;

    return sign * acos ( scalar );
}


float MainWindow::angleEE(MyMesh* _mesh, int vertexID,  int faceID)
{
    FaceHandle fh = _mesh->face_handle ( faceID );
    VertexHandle vh = _mesh->vertex_handle ( vertexID );
    std::vector<VertexHandle> vertexes;

    MyMesh::FaceVertexIter fh_v = _mesh->fv_iter(fh);

    for(; fh_v.is_valid(); ++fh_v) {
        VertexHandle current = *fh_v;
        if( current.idx() != vertexID )
            vertexes.push_back ( current );
    }

    float valuesAB[3];
    valuesAB[0] = _mesh->point(vertexes[0])[0] - _mesh->point(vh)[0];
    valuesAB[1] = _mesh->point(vertexes[0])[1] - _mesh->point(vh)[1];
    valuesAB[2] = _mesh->point(vertexes[0])[2] - _mesh->point(vh)[2];

    float valuesAC[3];
    valuesAC[0] = _mesh->point(vertexes[1])[0] - _mesh->point(vh)[0];
    valuesAC[1] = _mesh->point(vertexes[1])[1] - _mesh->point(vh)[1];
    valuesAC[2] = _mesh->point(vertexes[1])[2] - _mesh->point(vh)[2];

    VectorT<float, 3> normalizedAB(VectorT<float, 3>(valuesAB).normalize());
    VectorT<float, 3> normalizedAC(VectorT<float, 3>(valuesAC).normalize());

    return acos ( normalizedAB | normalizedAC );
}


void MainWindow::H_Curv(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++) {
        MyMesh::VertexHandle current = *curVert;
        float val = 0.0f;

        for (MyMesh::VertexEdgeIter currentEdge = _mesh->ve_iter ( current ); currentEdge.is_valid(); currentEdge++)
        {
            MyMesh::EdgeHandle eh = *currentEdge;
            MyMesh::HalfedgeHandle heh0 = _mesh->halfedge_handle(eh, 0);
            MyMesh::HalfedgeHandle heh1 = _mesh->halfedge_handle(eh, 1);

            FaceHandle fh0 = _mesh->face_handle(heh0);
            FaceHandle fh1 = _mesh->face_handle(heh1);

            // Si l'arête est en bordure, on ne traite qu'une face
            if ( fh1.idx ( ) > _mesh->n_faces ( ) )
                fh1 = fh0;

            // Détermine l'autre sommet
            int vertex2ID = _mesh->to_vertex_handle(heh1).idx();
            if (vertex2ID == current.idx ( ) )
                vertex2ID = _mesh->to_vertex_handle(heh0).idx();

            // ||e_ij||
            //VectorT<float,3> e = _mesh->point(_mesh->vertex_handle(vertex2ID)) - _mesh->point(vh);
            OpenMesh::Vec3f currentOppVector = _mesh->point ( _mesh->vertex_handle ( vertex2ID ) ) - _mesh->point ( current );

            OpenMesh::Vec3f normal0 ( _mesh->normal ( fh0 ) );
            OpenMesh::Vec3f normal1 ( _mesh->normal ( fh1 ) );

            if ( ( ( normal0 % normal1 ) | currentOppVector ) < 0 )
            {
                FaceHandle tempF = fh0;
                fh0 = fh1;
                fh1 = tempF;
            }

            val += currentOppVector.norm ( ) * angleFF ( _mesh , fh0.idx ( ) , fh1.idx ( ) , current.idx() , vertex2ID );
        }
        val /= ( 4 * baryArea ( _mesh , current.idx ( ) ) );
        _mesh->data ( current ).value = val;
    }
}

void MainWindow::K_Curv(MyMesh* _mesh)
{
    for ( MyMesh::VertexIter curVert = _mesh->vertices_begin() ; curVert!=_mesh->vertices_end() ; ++curVert ) {
        VertexHandle current = *curVert;

        float area = baryArea ( _mesh , current.idx() );
        float angleEESum = 0;

        MyMesh::VertexFaceIter vf = _mesh->vf_iter ( current );
        for ( ; vf.is_valid ( ) ; ++vf ) {
            FaceHandle currentFace = *vf;
            angleEESum += angleEE ( _mesh , current.idx ( ) , currentFace.idx ( ) );
        }

        _mesh->data ( current ).value = ( 1 / area ) * ( 2 * M_PI - angleEESum );
    }
}

//VectorT<unsigned, 37> MainWindow::classifyAngleFF (MyMesh* _mesh)
//{
//    VectorT<unsigned, 37> result;
//    for (unsigned i = 0 ; i < 37 ; ++i)
//        result[i] = 0;

//    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin() ; curEdge != _mesh->edges_end() ; ++curEdge)
//    {
//        MyMesh::EdgeHandle eh = *curEdge;
//        MyMesh::HalfedgeHandle heh0 = _mesh->halfedge_handle(eh, 0);
//        MyMesh::HalfedgeHandle heh1 = _mesh->halfedge_handle(eh, 1);

//        MyMesh::FaceHandle fh0 = _mesh->face_handle(heh0);
//        MyMesh::FaceHandle fh1 = _mesh->face_handle(heh1);

//        OpenMesh::Vec3f normal0 (_mesh->normal ( fh0 ) );
//        OpenMesh::Vec3f normal1 (_mesh->normal ( fh1 ) );
//        double angle = (normal0 | normal1);

//        if (angle > 1)
//            angle = 1;

//        qDebug() << unsigned((acos(angle)*180/M_PI)/10) << "?";
//        result[unsigned((acos(angle)*180/M_PI)/10)] ++;
//        qDebug() << unsigned((acos(angle)*180/M_PI)/10) << "!";

//    }

//    for (unsigned i = 0 ; i < 37 ; ++i)
//    {
//        qDebug() << "Angles à " << i*10 << "° :" << result[i];
//    }

//    return result;
//}

VectorT<float,2> MainWindow::minmaxAreaFace (MyMesh* _mesh)
{
    VectorT<float,2> minmax;
    minmax[0] = INT_MAX;
    minmax[1] = 0;
    for (int i = 0 ; i < _mesh->n_faces() ; ++i)
    {
        float area = faceArea(_mesh,i);
        if (area < minmax[0])
            minmax[0] = area;
        if (area > minmax[1])
            minmax[1] = area;
    }
    return minmax;
}

void MainWindow::displayFaceAreaFreq (MyMesh* _mesh)
{
    VectorT<float,2> minmax = minmaxAreaFace(_mesh);

    qDebug() << "min : " << minmax[0] << " ; max = " << minmax[1];

    VectorT<float,10> pallier;
    VectorT<int,10> cmpt;
    for (unsigned i = 0 ; i < 10 ; ++i)
    {
        float fac = 0.1f*(i+1);
        pallier[i] = minmax[0] + (minmax[1]-minmax[0])*fac;
        cmpt[i] = 0;
    }
    for (int i = 0 ; i < _mesh->n_faces() ; ++i)
    {
        bool found = false;
        float area = faceArea(_mesh,i);
        for (unsigned i = 0 ; i < 9 && !found ; ++i)
            if (pallier[i] < area  && area <= pallier[i+1])
            {
                found = true;
                cmpt[i+1] ++;
            }
        if (!found)
            cmpt[0]++;
    }

    qDebug() << "Pallier " << 0 << " : [" << minmax[0] << "," << pallier [0] << "] =" << cmpt[0];

    for (unsigned i = 0 ; i < 9 ; ++i)
        qDebug() << "Pallier " << i+1 << " : ]" << pallier[i] << "," << pallier [i+1] << "] =" << cmpt[i+1];
}

bool MainWindow::isFaceTriangle(MyMesh* _mesh, int faceID)
{
    FaceHandle fh = _mesh->face_handle (faceID);
    int count = 0;
    for (MyMesh::FaceVertexIter curVer = _mesh->fv_iter(fh) ; curVer.is_valid() ; curVer++)
        count++;
    return (count==3);
}

bool MainWindow::isAllFaceTriangle(MyMesh* _mesh)
{
    bool isTriangle = true;
    for (unsigned i = 0 ; i < _mesh->n_faces() && isTriangle ; i++)
        isTriangle = isFaceTriangle(_mesh, i);
    return isTriangle;
}

bool MainWindow::hasFaceNeighboors(MyMesh* _mesh, int faceID)
{
    FaceHandle fh = _mesh->face_handle (faceID);
    int count = 0;
    for (MyMesh::FaceFaceIter curFace = _mesh->ff_iter(fh) ; curFace.is_valid() && count == 0 ; curFace++)
        count++;
    return count>0;
}

bool MainWindow::areThereSingleFaces(MyMesh* _mesh)
{
    bool single = false;
    for (unsigned i = 0 ; i < _mesh->n_faces() && !single ; i++)
        single = !hasFaceNeighboors(_mesh, i);
    return single;
}

bool MainWindow::hasHalfEdgeAFace(MyMesh* _mesh, int halfEdgeID)
{
    HalfedgeHandle heh = _mesh->halfedge_handle(halfEdgeID);

}

VectorT <float,6> MainWindow::boundingBox3D(MyMesh* _mesh)
{
    VectorT <float,6> minmax;
    for(int i=0;i<3;i++){
        minmax[i]=_mesh->point(*_mesh->vertices_begin())[i];
        minmax[i+3]=_mesh->point(*_mesh->vertices_begin())[i];
    }
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        for(int i=0;i<3;i++){
            minmax[i]=std::min(minmax[i], _mesh->point(*curVert)[i]);
            minmax[i+3]=std::max(minmax[i+3], _mesh->point(*curVert)[i]);
        }
    }
    return minmax;
}


float MainWindow::totFaceArea(MyMesh* _mesh)
{
    float totFace = 0.0;
    for (unsigned int i = 0; i < _mesh->n_faces() ; i++)
    {
        totFace += faceArea(_mesh, i);
    }
    return totFace;
}

void calcul_barycentre(MyMesh* _mesh)
{
    float totX = 0.0;
    float totY = 0.0;
    float totZ = 0.0;
    float i = 0;
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        VertexHandle vh = *curVert;
        totX += _mesh->point(vh)[0];
        totY += _mesh->point(vh)[1];
        totZ += _mesh->point(vh)[2];
        i++;
    }
    totX /= i;
    totY /= i;
    totZ /= i;

    qDebug() << "Centre de gravité : (" << totX << ", " << totY << ", " << totZ << ")";
}


Vec3f normale_sommet (MyMesh* _mesh, int vertexID)
{
    return OpenMesh::Vec3f (_mesh->normal ( _mesh->vertex_handle(vertexID) ) );
}

Vec3f normale_face (MyMesh* _mesh, int faceID)
{
    return OpenMesh::Vec3f (_mesh->normal ( _mesh->face_handle(faceID) ) );
}

void afficher_normales_faces_sommets (MyMesh* _mesh)
{
    for (unsigned int i = 0; i < _mesh->n_faces() ; i++)
    {
        OpenMesh::Vec3f norface = normale_face(_mesh, i);
        qDebug() << "Normale face " << i << " : (" << norface[0] << ", " << norface[1] << ", " <<  norface[2] << ")";
    }

    for (unsigned int i = 0; i < _mesh->n_vertices() ; i++)
    {
        OpenMesh::Vec3f norvert = normale_sommet(_mesh, i);
        qDebug() << "Normale sommet " << i << " : (" << norvert[0] << ", " << norvert[1] << ", " <<  norvert[2] << ")";
    }
}




float MainWindow::ecart_angulaire_max(MyMesh* _mesh, int vertexID){
    VertexHandle v_it = _mesh->vertex_handle(vertexID);
    OpenMesh::Vec3f norvert = normale_sommet(_mesh, vertexID);
    OpenMesh::Vec3f norface;
    float maxi=0,temp;
    for (MyMesh::VertexFaceIter curFace = _mesh->vf_iter(v_it); curFace.is_valid(); curFace ++)
    {
           norface= normale_face(_mesh,(*curFace).idx());
           temp=acos(norface | norvert);
           maxi= std::max(maxi,temp);

    }
    return maxi;
}


void MainWindow::Deviation_normales(MyMesh *_mesh){
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        VertexHandle vh = *curVert;
        float value = ecart_angulaire_max(_mesh, vh.idx());
        _mesh->data(vh).value = value;
    }
}

/* **** fin de la partie à compléter **** */



/* **** début de la partie boutons et IHM **** */
void MainWindow::on_pushButton_H_clicked()
{
    H_Curv(&mesh);
    displayMesh(&mesh, true); // true permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

void MainWindow::on_pushButton_K_clicked()
{
    K_Curv(&mesh);
    calcul_barycentre(&mesh);
    afficher_normales_faces_sommets (&mesh);
    displayMesh(&mesh, true); // true permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

void MainWindow::on_pushButton_D_clicked()
{
    Deviation_normales(&mesh);
    displayMesh(&mesh, true); // true permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

/*
    Cette fonction est à utiliser UNIQUEMENT avec le fichier testAngleArea.obj
    Elle est appelée par le bouton "Test angles/aires"

    Elle permet de vérifier les fonctions faceArea, angleFF et angleEE.
    Elle doit afficher :

    Aire de la face 0 : 2
    Aire de la face 1 : 2
    Angle entre les faces 0 et 1 : 1.5708
    Angle entre les faces 1 et 0 : -1.5708
    Angle au sommet 1 sur la face 0 : 0.785398
*/

void MainWindow::on_pushButton_angleArea_clicked()
{
    qDebug() << "Aire de la face 0 :" << faceArea(&mesh, 0);
    qDebug() << "Aire de la face 1 :" << faceArea(&mesh, 1);

    qDebug() << "Angle entre les faces 0 et 1 :" << angleFF(&mesh, 0, 1, 1, 2);
    qDebug() << "Angle entre les faces 1 et 0 :" << angleFF(&mesh, 1, 0, 1, 2);

    qDebug() << "Angle au sommet 1 sur la face 0 :" << angleEE(&mesh, 1, 0);
    qDebug() << "Angle au sommet 3 sur la face 1 :" << angleEE(&mesh, 3, 1);
}

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);

    if (isAllFaceTriangle(&mesh))
        qDebug() << "Toutes les faces sont des triangles";
    else
        qDebug() << "Toutes les faces ne sont pas des triangles";

    if (areThereSingleFaces(&mesh))
        qDebug() << "Il y a au moins une face isolée.";
    else
        qDebug() << "Toutes les faces ont un voisin au moins.";

    qDebug() << "Nombre de faces : " << mesh.n_faces();
    qDebug() << "Nombre de sommets : " << mesh.n_vertices();

    VectorT <float,6> minmax;
    minmax=boundingBox3D(&mesh);

    qDebug() << "x_min: " << minmax[0] << "x_max: " << minmax[3];
    qDebug() << "y_min: " << minmax[1] << "y_max: " << minmax[4];
    qDebug() << "z_min: " << minmax[2] << "z_max: " << minmax[5];

    displayFaceAreaFreq(&mesh);

    classifyAngleFF(&mesh);

    qDebug() << "Aire totale du maillage : " << totFaceArea(&mesh);
}
/* **** fin de la partie boutons et IHM **** */

/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}
