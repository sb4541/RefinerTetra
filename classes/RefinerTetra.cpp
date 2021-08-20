#include "RefinerTetra.hpp"
#include "Output.hpp"
#include "Eigen"
#include "Test.h"
#include "algorithm"

using namespace MainApplication;

namespace GeDiM
{
    RefinerTetra::RefinerTetra()
    {
        meshPointer = NULL;
        idCellToRefine = vector<unsigned int>{};
        idEdgesCut = list<unsigned int>{};
        //check_presence = vector<bool>{};
        SetCounter();
    }



//****************************************************************************************************
//****************************************************************************************************
    RefinerTetra::~RefinerTetra()
    {
        meshPointer = NULL;
        idEdgesCut.clear();
        idCellToRefine.clear();
        //check_presence.clear();
    }





//****************************************************************************************************
//****************************************************************************************************
    /*void RefinerTetra::SetVector()
    {
        unsigned int numberofEdges = meshPointer->NumberOfEdges();
        check_presence.reserve(numberofEdges);
        for(int position = 0; position < numberofEdges; ++position)
            check_presence.push_back(false);
    }*/


//****************************************************************************************************
//****************************************************************************************************
    const Output::ExitCodes RefinerTetra::FindQualityCell(GenericCell &cell, double &quality, bool &check_problem)
    {
        //Calcola l'aspect ratio come rapporto tra lato più lungo e minima altezza

        /*Vector3d centroid;
        for(unsigned int position = 0; position < cell.NumberOfFaces(); ++position)
        {
            unsigned int idFace = cell.Face(position)->Id();
            GenericFace& face = *meshPointer->Face(idFace);
            face.ComputeGeometricalProperties();
        }
        cell.Centroid(centroid);
        double tolerance = 1.0E-12;
        //'bool' per segnalare un problema descritto alla fine del metodo
        check_problem = false;

        //Si calcola il raggio esterno
        double external_radius = 0.0;
        cell.Radius(external_radius, centroid);*/
        //Se il raggio esterno è troppo piccolo, il metodo precedente genera un errore e non si prosegue (il raggio interno sarà ancora più piccolo)

        //RecoverConformity();
        double tolerance = 1.0E-12;

        //unsigned int idLongestEdge;
        unsigned int NumberOfEdges = cell.NumberOfEdges();
        double current_length = 0.0;
        for(int edg = 0; edg < NumberOfEdges; ++edg)
        {
            unsigned int idTemp = cell.Edge(edg)->Id();
            GenericEdge& edge = *meshPointer->Edge(idTemp);
            double length = (edge.Point(0)->Coordinates() - edge.Point(1)->Coordinates()).squaredNorm();
            if(current_length < length)
            {
                //idLongestEdge = idTemp;
                current_length = length;
            }
        }
        //FindLongestEdge(cell, idLongestEdge);
        //GenericEdge& LongestEdge = *meshPointer->Edge(idLongestEdge);
        double numerator = sqrt(current_length);

        //Si calcola il raggio interno con Intersector2D1D (retta da centroide perpendicolare ad ogni faccia; lo si setta come il raggio esterno perché questo è più grande di quello interno)
        double denominator = numerator + 1.0;
        map<unsigned int, unsigned int> points;
        for(int position = 0; position < cell.NumberOfPoints(); ++position)
        {
            unsigned int idTemp = cell.Point(position)->Id();
            points.insert(pair<unsigned int, unsigned int>(idTemp, 0));
        }
        for (int position = 0; position < cell.NumberOfFaces(); ++position)
        {
            GenericFace &face = *meshPointer->Face(cell.Face(position)->Id());
            face.ComputeNormal();
            const Vector3d normal = face.Normal();
            double planeTranslation = normal.dot(face.Point(0)->Coordinates());
            Intersector2D1D intersector;
            for(int point = 0; point < face.NumberOfPoints(); ++point)
            {
                unsigned int idTemp = face.Point(point)->Id();
                map<unsigned int, unsigned int>::iterator it = points.find(idTemp);
                it->second = 1;
            }
            unsigned int idOpposite;
            for(int point = 0; point < cell.NumberOfPoints(); ++point)
            {
                unsigned int idTemp = cell.Point(point)->Id();
                map<unsigned int, unsigned int>::iterator it = points.find(idTemp);
                if(it->second == 0)
                {
                    idOpposite = it->first;
                    break;
                }
            }
            GenericPoint& OppositePoint = *meshPointer->Point(idOpposite);
            intersector.SetLine(OppositePoint.Coordinates(), normal);
            intersector.SetPlane(normal, planeTranslation);
            intersector.ComputeIntersection();
            const Vector3d intersectionPoint = intersector.IntersectionPoint();
            double squaredDistance = (OppositePoint.Coordinates() - intersectionPoint).squaredNorm();
            if (denominator > squaredDistance)
                denominator = squaredDistance;
            map<unsigned int, unsigned int>::iterator it = points.begin();
            for(; it != points.end(); ++it)
            {
                it->second = 0;
            }
        }
        denominator = sqrt(denominator);
        //Se il raggio interno è troppo piccolo (ma non quello grande, se no metodo si ferma prima), si decide di settarlo a zero e indicare questa situazione con un 'bool' (vedere l'uso di esso in 'FindWorstQualityCells')
        if (fabs(denominator) < tolerance)
        {
            denominator = 0.0;
            check_problem = true;
            return Output::Success;
        }

        quality = numerator / denominator;

        return Output::Success;
    }



//****************************************************************************************************
//****************************************************************************************************
    inline const Output::ExitCodes RefinerTetra::FindQualityCell_2(GenericCell &cell, double &quality)
    {
        quality = 0.0;
        //Si è scelto il rapporto tra massimo e minimo lato come valutazione di qualità
        //(non è una valutazione generale, ma è usato soltanto in 'FindWorstQualityCells' per ordinare le celle in base ad esso)
        quality = cell.RatioMaxMinEdge();

        return Output::Success;
    }



//****************************************************************************************************
//****************************************************************************************************
    const Output::ExitCodes RefinerTetra::FindWorstQualityCells(multimap<double, unsigned int> &WorstQualityCells, unsigned int &num_Cells, double& mean)
    {
        //Prima ci si assicura che ci sia conformità, perché questo metodo viene chiamato su mesh conforme
        RecoverConformity();
        multimap<double, unsigned int> Current_WorstQualityCells;
        unsigned int numberofCells = meshPointer->NumberOfCells();
        unsigned int divisor = 0;
        mean = 0.0;

        for (int id_cell = 0; id_cell < numberofCells; ++id_cell)
        {
            GenericCell &cell_toStudy = *meshPointer->Cell(id_cell);
            //Guardi solo le celle attive
            if (!cell_toStudy.IsActive())
            {
                continue;
            }
            else
            {
                double quality = 0.0;
                bool check = false;
                FindQualityCell(cell_toStudy, quality, check);
                if(check == false)
                {
                    ++divisor;
                    mean += quality;
                    //Rese negative perché così quelle con qualità pessima (dato da alto rapporto MaxMinEdge) compaiono per prime visto che l'ordine in mappa è decrescente
                    quality = -quality;
                    Current_WorstQualityCells.insert(std::pair<double, unsigned int>(quality, id_cell));
                }

 /*               FindQualityCell_2(cell_toStudy, quality);
                ++divisor;
                mean += quality;
                //Rese negative perché così quelle con qualità pessima (dato da alto rapporto MaxMinEdge) compaiono per prime visto che l'ordine in mappa è decrescente
                quality = -quality;
                Current_WorstQualityCells.insert(std::pair<double, unsigned int>(quality, id_cell));
*/
//
            }

        }
        mean = (double) mean/divisor;

        //Si vuole soltanto un certo numero, ossia 'num_Cells' di celle, quindi copio soltanto le prime .....
        std::multimap<double, unsigned int>::iterator it = Current_WorstQualityCells.begin();
        //Metto entrambe condizioni in 'for' perché: in 'PrintStatistics' viene inserito di default l'intera foresta di celle, mentre talvolta
        //si desidera le prime ... (numero dato da 'num_Cells')
        for (int id_element = 0; it != Current_WorstQualityCells.end() && id_element < num_Cells; ++id_element, ++it)
        {
            WorstQualityCells.insert(std::pair<double, unsigned int>(it->first, it->second));
        }

        //OPPURE
        //si usa 'FindQualityCell'
        /*
        RecoverConformity();
        multimap<double, unsigned int> Current_WorstQualityCells;
        //'degenerates' per gestire le celle che hanno raggio interno troppo piccolo (le si metterà all'inizio della mappa di celle supponendo siano le peggiori)
        vector<unsigned int> degenerates;

        for(int id_cell = 0; id_cell < meshPointer->NumberOfCells(); ++id_cell)
        {
            GenericCell& cell_toStudy = *meshPointer->Cell(id_cell);
            //Guardi solo le celle attive
            if(!cell_toStudy.IsActive())
            {
                continue;
            }
            else
            {
               bool check_problem = false;
               double quality = 0.0;
               FindQualityCell(cell_toStudy, quality, check_problem);
               if(check_problem)
               {
                    degenerates.push_back(id_cell);
               }
               else
               {
                   //Rese negative perché così quelle con qualità pessima (dato da alto rapporto MaxMinEdge) compaiono per prime visto che l'ordine in mappa è decrescente
                   quality = -1.0 * quality;

                   Current_WorstQualityCells.insert(std::pair<double, unsigned int>(quality, id_cell));
               }
            }

        }

        //Si setta la qualità peggiore per i 'degenerates' così da metterli per primi (in linea con la strategia descritta precedentemente)
        double quality_default = Current_WorstQualityCells.begin()->first - 1;
        for(int position = 0; position < degenerates.size(); ++position)
        {
            Current_WorstQualityCells.insert(std::pair<double, unsigned int>(quality_default, degenerates[position]));
        }

        std::multimap<double, unsigned int>::iterator it = Current_WorstQualityCells.begin();
        for(int id_element = 0; id_element < num_Cells; ++id_element, ++it)
        {
            WorstQualityCells.insert(std::pair<double, unsigned int>(it->first, it->second));
        }
        */


        return Output::Success;
    }


//****************************************************************************************************
//****************************************************************************************************
    const Output::ExitCodes RefinerTetra::PrintStatistics(unsigned int NumberOfCells)
    {
        RecoverConformity();
        printf("%d cuts have been carried out\n", refiner_counter);

        multimap<double, unsigned int> WorstQualityCells;
        unsigned int num_Cells = meshPointer->NumberOfCells();
        double mean = 0.0;
        FindWorstQualityCells(WorstQualityCells, NumberOfCells, mean);

        std::multimap<double, unsigned int>::iterator it = WorstQualityCells.begin();
        printf("Worst quality is from %d: %lf\n\n", it->second, (-1.0)*it->first);

        printf("These are the %d worst cells:\n", NumberOfCells);
        for(; it != WorstQualityCells.end(); ++it)
        {
            printf("%d: %lf\n", it->second, (-1.0)*it->first);
        }
        printf("\n");

        printf("The mean quality is: %lf\n\n", mean);

        return Output::Success;
    }




//****************************************************************************************************
//****************************************************************************************************
    const Output::ExitCodes RefinerTetra::FindMaxFace(GenericCell &cell, unsigned int &idMaxFace)
    {
        //idMaxFace conterrà l'id globale della faccia di cell con area maggiore
        double MaxAreaFace = 0.0;
        idMaxFace = meshPointer->NumberOfFaces();
        int number_of_faces = cell.NumberOfFaces();  //da usare in ciclo for seguente

        for (int position = 0; position < number_of_faces; ++position)
        {
            GenericFace &faceToTreat = *meshPointer->Face(cell.Face(position)->Id());
            faceToTreat.ComputeGeometricalProperties();         //viene calcolato il membro measure che è l'area della faccia

            if (faceToTreat.Measure() > MaxAreaFace)
            {
                idMaxFace = position;
                MaxAreaFace = faceToTreat.Measure();
            }
        }
        idMaxFace = cell.Face(idMaxFace)->Id();

        if (idMaxFace == meshPointer->NumberOfFaces())
        {
            Output::PrintErrorMessage("Biggest face was not found\n", false);
            return Output::GenericError;
        }

#if VERBOSITY == UNIT_TEST
        UnitTest::Test_MaxFace(cell, idMaxFace);
#endif
        return Output::Success;
    }


//****************************************************************************************************
//****************************************************************************************************
    const Output::ExitCodes RefinerTetra::FindLongestEdge(const GenericCell &cell, unsigned int &idLongestEdge)
    {
        //Avrà id globale di lato più lungo di intera cella
        idLongestEdge = meshPointer->NumberOfEdges();
        unsigned int idFace = meshPointer->NumberOfFaces();
        double length_longest_edge = 0.0;
        int number_of_faces = cell.NumberOfFaces();

        for (int position = 0; position < number_of_faces; ++position)
        {
            unsigned int idTemp = cell.Face(position)->Id();
            GenericFace& face = *meshPointer->Face(idTemp);
            int number_of_edges = face.NumberOfEdges();
            for(int edge = 0; edge < number_of_edges; ++edge)
            {
                unsigned int _idTemp = face.Edge(edge)->Id();
                GenericEdge& edgeToTreat = *meshPointer->Edge(_idTemp);
                //id in Point sono 0 e 1 come fatto a riga 2573 in GenericMesh.cpp
                double current_value = (edgeToTreat.Point(0)->Coordinates() - edgeToTreat.Point(1)->Coordinates()).squaredNorm();
                if (current_value > length_longest_edge)
                {
                    idLongestEdge = edge;
                    idFace = position;
                    length_longest_edge = current_value;
                }
                //Nel caso in cui ci siano due lati di uguale lunghezza, prendi quello che è già stato tagliato se c'è, per non avere errore in 'CutTetra'
                else if (current_value == length_longest_edge && !edgeToTreat.IsActive())
                {
                    //Controllo per capire se questo è il lato che è stato tagliato (attraverso l'id del punto medio di taglio)
                    //Può succedere che una faccia sia tagliata lungo un lato ma anche un altro lato non sia più attivo (io non voglio quest'ultimo)
                    if(!face.IsActive())
                    {
                        unsigned int idChild = edgeToTreat.Child(0)->Id();
                        unsigned int idMiddlePoint = meshPointer->Edge(idChild)->Point(0)->Id();
                        idChild = face.Child(0)->Id();
                        unsigned int idMiddlePoint_toCompare = meshPointer->Face(idChild)->Point(0)->Id();
                        if(idMiddlePoint == idMiddlePoint_toCompare)
                        {
                            idLongestEdge = edge;
                            idFace = position;
                            //length_longest_edge = current_value;
                        }
                    }
                    else if(face.IsActive() && cell.Face(idFace)->IsActive() && cell.Face(idFace)->Edge(idLongestEdge)->IsActive() && !edgeToTreat.IsActive())
                    {
                        idLongestEdge = edge;
                        idFace = position;
                        //length_longest_edge = current_value;
                    }
                }
            }
        }
        idLongestEdge = cell.Face(idFace)->Edge(idLongestEdge)->Id();

        if (idLongestEdge == meshPointer->NumberOfEdges())
        {
            Output::PrintErrorMessage("Longest edge was not found\n", false);
            return Output::GenericError;
        }

        return Output::Success;
    }


//****************************************************************************************************
//****************************************************************************************************

    //SOLO LE FACCE E I TETRAEDRI POSSONO AVERE 'NULL' NEI VETTORI DEI TETRAEDRI SE SONO DI BORDO (GLI ALTRI OGGETTI, NON HANNO 'NULL')
    //SI E' DECISO DI SOSTITUIRE I FIGLI CON IL PADRE (DI QUALSIASI OGGETTO) NEI VICINI DEGLI OGGETTI DA CUI SONO FORMATI I FIGLI TRANNE CHE PER IL LATO PII' LUNGO
    //CHE SARA' USATO PER LA CONFORMITA'

    const Output::ExitCodes RefinerTetra::CutTetra(GenericCell &cell, bool& to_add)
    {
        ++refiner_counter;

        //L'idea è di prendere il lato più lungo e tagliare la cella attraverso il punto
        //punto medio del lato più lungo e i due punti che non appartengono al lato maggiore (id_first_opposite_point, id_second_opposite_point)
        unsigned int id_first_opposite_point = meshPointer->NumberOfPoints();
        unsigned int id_second_opposite_point = meshPointer->NumberOfPoints();

        unsigned int idLongestEdge = 0;
        FindLongestEdge(cell, idLongestEdge);
        GenericEdge &LongestEdge = *meshPointer->Edge(idLongestEdge);
        unsigned int id_first_point = LongestEdge.Point(0)->Id();
        unsigned int id_second_point = LongestEdge.Point(1)->Id();
        GenericPoint &first_point = *meshPointer->Point(id_first_point);
        GenericPoint &second_point = *meshPointer->Point(id_second_point);

        int number_vertices_cell = cell.NumberOfPoints();
        unsigned int flag = 0;   //per sapere se hai già settato id_first_opposite_point e per uscire poi dal ciclo for (vedi ciclo for)
        //ciclo for per trovare i due vertici che non sono estremi della faccia maggiore
        for (int position = 0; position < number_vertices_cell; ++position)
        {
            if (flag == 2)
            {
                break;
            }
            if (cell.Point(position)->Id() != id_first_point && cell.Point(position)->Id() != id_second_point && flag == 0)
            {
                id_first_opposite_point = position;
                ++flag;
            }
            else if (cell.Point(position)->Id() != id_first_point && cell.Point(position)->Id() != id_second_point && flag == 1)
            {
                id_second_opposite_point = position;
                ++flag;
            }
        }
        id_first_opposite_point = cell.Point(id_first_opposite_point)->Id();
        id_second_opposite_point = cell.Point(id_second_opposite_point)->Id();

        if (flag != 2)
        {
            Output::PrintErrorMessage("Opposite vertices were not found\n", false);
            return Output::GenericError;
        }

        GenericPoint &first_opposite_point = *meshPointer->Point(id_first_opposite_point);
        GenericPoint &second_opposite_point = *meshPointer->Point(id_second_opposite_point);

//****************************************************************************************************
        unsigned int id_first_face = meshPointer->NumberOfFaces();         //conterrà first_opposite_point
        unsigned int id_second_face = meshPointer->NumberOfFaces();        //conterrà second_opposite_point
        //ciclo costruito in modo che le due facce che voglio (ossia quelle da tagliare, ossia quelle con il lato più lungo e rispettivamente 'first_opposite_point'
        // e 'second_opposite_point')
        //avranno counter = 4 e counter = 5, mentre le altre due avranno counter = 6
        flag = 0;
        for (int position = 0; position < cell.NumberOfFaces(); ++position)
        {
            if(flag == 2)
                break;
            unsigned int counter = 0;
            GenericFace &face = *meshPointer->Face(cell.Face(position)->Id());

            //NOTA NOTA NOTA
            unsigned int r = 0;

            for (int point = 0; point < face.NumberOfPoints(); point++)
            {
                unsigned int face_point_id = face.Point(point)->Id();
                r = face_point_id;
                if (face_point_id == id_first_point || face.Point(point)->Id() == id_second_point)
                {
                    ++counter;
                }
                if (face_point_id == id_first_opposite_point)
                    counter = counter + 2;
                if (face_point_id == id_second_opposite_point)
                    counter = counter + 3;
            }
            if (counter == 4)
            {
                ++flag;
                id_first_face = face.Id();
            }
            if (counter == 5)
            {
                id_second_face = face.Id();
                ++flag;
            }
            if(counter != 4 && counter != 5 && counter != 6)
            {
                Output::PrintErrorMessage("Face with id %d in cell with id %d has a wrong side\n", false, r, face.Id(), cell.Id());
                return Output::GenericError;
            }

        }

        if (flag != 2)
        {
            Output::PrintErrorMessage("Faces of id %d with the longest edge were not found\n", false, cell.Id());
            return Output::GenericError;
        }

        flag = 0;
        GenericFace &first_face = *meshPointer->Face(id_first_face);
        GenericFace &second_face = *meshPointer->Face(id_second_face);

//****************************************************************************************************
        //Si creano lati figli del lato più lungo (se non tagliato precedentemente) e si crea nuovo punto e si settano i suoi membri, se no li si prendono
        //Se è attivo, allora nessuna delle due facce di cui fa parte è stata tagliata (serve perché sai che il nuovo punto apparterrà soltanto ai due tetraedri figli
        // per ora, che non devi recuperare, ma devi creare ora; inoltre, se il lato è attivo, allora sai che devi proprio fare il taglio completo e non andare a prendere
        //facce figlie)

        if (LongestEdge.IsActive())
        {
            to_add = true;

            //si trova punto medio di lato più lungo e lo si inserisce nella mesh
            Vector3d halfpoint_longest_edge;
            halfpoint_longest_edge = 0.5 * (first_point.Coordinates() - second_point.Coordinates());
            halfpoint_longest_edge = second_point.Coordinates() + halfpoint_longest_edge;
            GenericPoint &newPoint = *meshPointer->CreatePoint();
            meshPointer->AddPoint(&newPoint);

            //Si riempiono membri del nuovo punto
            newPoint.SetCoordinates(static_cast<const Vector3d> (halfpoint_longest_edge));
            //'+2' perché ci saranno le due facce figlie e i due tetraedri figli (apparterrà indicativamente a stesso numero di facce e tetraedri: saranno tutte le
            //facce figli date da taglio di questo tetraedro e quelli con lato più lungo in comune)
            //LA PRESENZA DI OGGETTI NON PIU' ATTIVI IN VETTORE DI VICINI E' GESTITO IN CONFORMITA'
            newPoint.InitializeFaces(LongestEdge.NumberOfFaces() + 2);
            newPoint.InitializeCells(LongestEdge.NumberOfCells() + 2);


            //Si creano lati figli del lato più lungo
            LongestEdge.SetState(false);
            LongestEdge.InitializeChilds(2);
            GenericEdge &newEdge_1 = *(meshPointer->CreateEdge());
            meshPointer->AddEdge(&newEdge_1);
            newEdge_1.AddPoint(&newPoint);
            newEdge_1.AddPoint(&first_point);
            LongestEdge.AddChild(&newEdge_1);
            newEdge_1.SetFather(&LongestEdge);
            newEdge_1.InheritPropertiesByFather();
            GenericEdge &newEdge_2 = *(meshPointer->CreateEdge());
            meshPointer->AddEdge(&newEdge_2);
            newEdge_2.AddPoint(&newPoint);
            newEdge_2.AddPoint(&second_point);
            LongestEdge.AddChild(&newEdge_2);
            newEdge_2.SetFather(&LongestEdge);
            newEdge_2.InheritPropertiesByFather();
            //check_presence.push_back(false);
            //check_presence.push_back(false);

            //Si inseriscono nuovi lati in nuovo punto (due aggiunti ora, due dopo: quelli creati in taglio di due facce)
            newPoint.InitializeEdges(2 + 2);
            newPoint.AddEdge(&newEdge_1);
            newPoint.AddEdge(&newEdge_2);

            //Si inseriscono nuovi lati in punti già esistenti (estremi di nuovi lati insieme al nuovo punto), non al posto di lato che si sta tagliando
            //perché appartengono ancora a quel lato, ma quel lato è semplicemente non attivo
            //No controllo se il puntatore è NULL, perché non ce ne sono nel vettore dei vicini (non in partenza, non si mettono con questo algoritmo)
            first_point.AddEdge(&newEdge_1);
            second_point.AddEdge(&newEdge_2);

            //Si inseriscono le celle a cui appartengono i nuovi lati (no '+2' perché ciascun lato sarà parte di un nuovo tetraedro),
            //ma non inseriti, sarà fatto per tetraedri figli dopo e per tetraedri vicini durante recupero di conformità, che richiama questo metodo)
            newEdge_1.InitializeCells(LongestEdge.NumberOfCells());
            newEdge_2.InitializeCells(LongestEdge.NumberOfCells());

            //Si vorrebbero inserire le facce nei nuovi lati (quattro dopo: quelle create in taglio di due facce; no '+2' come per celle sopra)
            //Ma le facce di altri tetraedri saranno poi tagliate e quindi saranno aggiunti i figli
            newEdge_1.InitializeFaces(LongestEdge.NumberOfFaces());
            newEdge_2.InitializeFaces(LongestEdge.NumberOfFaces());

            //I NUOVI PUNTI E LATI VENGONO AGGIUNTI IN NUOVI TETRAEDRI, E LORO SOTTO-COMPONENTI, IN QUESTA FUNZIONE (MENTRE IN ALTRI TETRAEDRI, E LORO SOTTO-COMPONENTI,
            //IN ALTRE PARTI DEL CODICE, QUANDO QUESTI SARANNO RESI CONFORMI)

            //****************************************************************************************************
            //Da qui, si inizia a tagliare first_face (so già, che se il lato più lungo non è mai stato tagliato (sono ancora dentro all''if')
            // allora neanche le facce lo sono; si devono riempire tutti i membri dei nuovi oggetti( per le due nuove facce: punti, lati e tretraedri; per nuovo lato:
            //punti, facce, tetraedri)
            if (!first_face.IsActive())
            {
                Output::PrintErrorMessage("Error in cutting face with global id %d of cell with global id %d: this face should be active\n",id_first_face, cell.Id(), false);
                return Output::GenericError;
            }

            first_face.SetState(false);
            first_face.InitializeChilds(2);

            //Si creano facce figlie della faccia che si sta tagliando
            GenericFace &first_newFace_1 = *(meshPointer->CreateFace());
            meshPointer->AddFace(&first_newFace_1);
            first_face.AddChild(&first_newFace_1);
            first_newFace_1.SetFather(&first_face);
            first_newFace_1.InheritPropertiesByFather();
            first_newFace_1.AllocateCells(2);
            GenericFace &first_newFace_2 = *(meshPointer->CreateFace());
            meshPointer->AddFace(&first_newFace_2);
            first_face.AddChild(&first_newFace_2);
            first_newFace_2.SetFather(&first_face);
            first_newFace_2.InheritPropertiesByFather();
            first_newFace_2.AllocateCells(2);
            //Le celle che condividono una faccia sono solo due (se quella vicina è stata tagliata, nota: lungo un altro lato su un'altra sua faccia,
            //non è un problema), ma non le (qui solo il tetraedro figlio di questa cella)
            //inseriamo perché non fatto qui, fatto poi nel ricreare la conformità per altro tetraedro (infatti
            //magari figlio di altro tetraedro avrà figli e nipoti e poi è più difficile risale a quale di essi ha questa faccia)

            //si aggiungono i lati alle nuove facce e le nuove facce ai lati
            first_newFace_1.InitializeEdges(first_face.NumberOfEdges());
            first_newFace_2.InitializeEdges(first_face.NumberOfEdges());
            first_newFace_1.AddEdge(&newEdge_1);
            first_newFace_2.AddEdge(&newEdge_2);
            for (unsigned int i = 0; i < first_face.NumberOfEdges(); i++)
            {
                unsigned int idTemp = first_face.Edge(i)->Id();
                //Controllo su lato per capire se è quello che si sta tagliando (ossia il più lungo)
                //no controllo su lato se attivo (perché anche se uno può non esserlo, perché tagliato in un tetraedro che condivide solo un lato,
                //non è di interesse qui)
                if (idTemp != idLongestEdge)
                {
                    GenericEdge &edge = *meshPointer->Edge(idTemp);
                    //Controllo su lato per capire a quale faccia deve appartenere (tenendo conto che 'first_newFace_1' avrà 'first_point', ......
                    if (edge.Point(0)->Id() == id_first_point || edge.Point(1)->Id() == id_first_point)
                    {
                        edge.AddFace(&first_newFace_1);
                        first_newFace_1.AddEdge(&edge);
                    }
                    else if (edge.Point(0)->Id() == id_second_point || edge.Point(1)->Id() == id_second_point)
                    {
                        edge.AddFace(&first_newFace_2);
                        first_newFace_2.AddEdge(&edge);
                    }
                    else
                        return Output::GenericError;
                }
            }
            //Prima non si era inserito la faccia massima perché sarebbe stata tagliata
            newEdge_1.AddFace(&first_newFace_1);
            newEdge_2.AddFace(&first_newFace_2);

            //Si aggiungono i punti alle facce e viceversa
            first_newFace_1.InitializePoints(3);
            first_newFace_2.InitializePoints(3);
            first_newFace_1.AddPoint(&newPoint);
            first_newFace_1.AddPoint(&first_point);
            //Per come strutturato, so che il terzo vertice (comune ai due figli) di questa faccia ('first_face') è 'first_opposite_point'
            first_newFace_1.AddPoint(&first_opposite_point);
            first_newFace_2.AddPoint(&newPoint);
            first_newFace_2.AddPoint(&second_point);
            first_newFace_2.AddPoint(&first_opposite_point);
            newPoint.AddFace(&first_newFace_1);
            newPoint.AddFace(&first_newFace_2);
            first_point.AddFace(&first_newFace_1);
            second_point.AddFace(&first_newFace_2);
            first_opposite_point.AddFace(&first_newFace_1);
            first_opposite_point.AddFace(&first_newFace_2);


            //Si crea nuovo lato che deve tagliare la faccia 'first_face'
            GenericEdge &first_newEdge = *(meshPointer->CreateEdge());
            meshPointer->AddEdge(&first_newEdge);
            //Aggiungiamo punto medio al nuovo lato e l'altro vertice (so che è 'first_opposite_point') e viceversa
            first_newEdge.AddPoint(&newPoint);
            newPoint.AddEdge(&first_newEdge);
            first_newEdge.AddPoint(&first_opposite_point);
            first_opposite_point.AddEdge(&first_newEdge);

            //check_presence.push_back(false);

            //Si aggiungono le nuove facce al nuovo lato (e viceversa: non si era ancora inserito il nuovo lato perché non ancora creato quando si trattavano
            //le facce figlie)
            first_newEdge.AddFace(&first_newFace_1);
            first_newEdge.AddFace(&first_newFace_2);
            first_newFace_1.AddEdge(&first_newEdge);
            first_newFace_2.AddEdge(&first_newEdge);

            //****************************************************************************************************
            //Da qui, si inizia a tagliare la seconda faccia che condivide il lato più lungo
            //Si sa già qual è: 'second_face'

            if (!second_face.IsActive())
            {
                //Non può essere non attiva, se no lato più lungo tagliato e quindi non attivo
                Output::PrintErrorMessage("Error in cutting face with global id %d of cell with global id %d: this face should be active\n",id_second_face, cell.Id(), false);
                return Output::GenericError;
            }

            second_face.SetState(false);
            second_face.InitializeChilds(2);

            //Si creano facce figlie
            GenericFace &second_newFace_1 = *(meshPointer->CreateFace());
            meshPointer->AddFace(&second_newFace_1);
            second_face.AddChild(&second_newFace_1);
            second_newFace_1.SetFather(&second_face);
            second_newFace_1.InheritPropertiesByFather();
            second_newFace_1.AllocateCells(2);
            //Le celle che condividono una faccia sono solo due ma non si inseriscono ora (il tetraedro figlio di questa cella sarà inserito dopo e il vicino
            //nel recupero della conformità: come per 'first_face')

            GenericFace &second_newFace_2 = *(meshPointer->CreateFace());
            meshPointer->AddFace(&second_newFace_2);
            second_face.AddChild(&second_newFace_2);
            second_newFace_2.SetFather(&second_face);
            second_newFace_2.InheritPropertiesByFather();
            second_newFace_2.AllocateCells(2);

            //si aggiungono i lati alle facce e le nuove facce ai lati
            second_newFace_1.InitializeEdges(3);
            second_newFace_2.InitializeEdges(3);
            second_newFace_1.AddEdge(&newEdge_1);
            second_newFace_2.AddEdge(&newEdge_2);
            for (unsigned int i = 0; i < second_face.NumberOfEdges(); i++)
            {
                unsigned int idTemp = second_face.Edge(i)->Id();
                //Controllo su lato per capire se è quello che si sta tagliando (ossia il più lungo)
                //no controllo su lato se attivo (perché anche se uno può non esserlo, perché tagliato in un tetraedro che condivide solo un lato,
                //non è di interesse qui, ma nel recupero della conformità)
                if (idTemp != idLongestEdge)
                {
                    GenericEdge &edge = *meshPointer->Edge(idTemp);
                    //Controllo su lato per capire a quale faccia deve appartenere (tenendo conto che 'first_newFace_1' avrà 'first_point', ......
                    if (edge.Point(0)->Id() == id_first_point || edge.Point(1)->Id() == id_first_point)
                    {
                        edge.AddFace(&second_newFace_1);
                        second_newFace_1.AddEdge(&edge);
                    }
                    else if (edge.Point(0)->Id() == id_second_point || edge.Point(1)->Id() == id_second_point)
                    {
                        edge.AddFace(&second_newFace_2);
                        second_newFace_2.AddEdge(&edge);
                    }
                    else
                        return Output::GenericError;
                }
            }
            newEdge_1.AddFace(&second_newFace_1);
            newEdge_2.AddFace(&second_newFace_2);

            //Si aggiungono i punti alle facce e viceversa
            second_newFace_1.InitializePoints(3);
            second_newFace_2.InitializePoints(3);
            second_newFace_1.AddPoint(&newPoint);
            second_newFace_1.AddPoint(&first_point);
            //Per come strutturato, so che il terzo vertice (comune ai due figli) di questa faccia ('second_face') è 'second_opposite_point'
            second_newFace_1.AddPoint(&second_opposite_point);
            second_newFace_2.AddPoint(&newPoint);
            second_newFace_2.AddPoint(&second_point);
            second_newFace_2.AddPoint(&second_opposite_point);
            newPoint.AddFace(&second_newFace_1);
            newPoint.AddFace(&second_newFace_2);
            first_point.AddFace(&second_newFace_1);
            second_point.AddFace(&second_newFace_2);
            second_opposite_point.AddFace(&second_newFace_1);
            second_opposite_point.AddFace(&second_newFace_2);


            //Si crea nuovo lato che deve tagliare la faccia 'second_face'
            GenericEdge &second_newEdge = *(meshPointer->CreateEdge());
            meshPointer->AddEdge(&second_newEdge);
            //Aggiungiamo punto medio al nuovo lato e l'altro vertice (so che è 'first_opposite_point') e viceversa
            second_newEdge.AddPoint(&newPoint);
            newPoint.AddEdge(&second_newEdge);
            second_newEdge.AddPoint(&second_opposite_point);
            second_opposite_point.AddEdge(&second_newEdge);

            //check_presence.push_back(false);

            //Si aggiungono le nuove facce al nuovo lato (e viceversa: non si era ancora inserito il nuovo lato perché non ancora creato quando si trattavano
            //le facce figlie)
            second_newEdge.AddFace(&second_newFace_1);
            second_newEdge.AddFace(&second_newFace_2);
            second_newFace_1.AddEdge(&second_newEdge);
            second_newFace_2.AddEdge(&second_newEdge);

            //****************************************************************************************************
            //Da qui si crea la terza nuova faccia

            GenericFace &third_newFace = *(meshPointer->CreateFace());
            meshPointer->AddFace(&third_newFace);
            third_newFace.AllocateCells(2);
            //Saranno inseriti i tetredri figli

            //Si aggiungono i punti alla faccia (si conoscono già) e viceversa
            third_newFace.InitializePoints(3);
            third_newFace.AddPoint(&newPoint);
            third_newFace.AddPoint(&first_opposite_point);
            third_newFace.AddPoint(&second_opposite_point);
            first_opposite_point.AddFace(&third_newFace);
            second_opposite_point.AddFace(&third_newFace);
            newPoint.AddFace(&third_newFace);

            //Si aggiungono i lati alla terza faccia e viceversa
            third_newFace.InitializeEdges(3);
            third_newFace.AddEdge(&first_newEdge);
            third_newFace.AddEdge(&second_newEdge);
            first_newEdge.AddFace(&third_newFace);
            second_newEdge.AddFace(&third_newFace);
            //Si deve andare a prendere il terzo lato che non è mai stato trovato (quello con 'first_opposite_point' e 'second_opposite_point'
            //come vertici)
            for (int position = 0; position < cell.NumberOfEdges(); ++position)
            {
                unsigned int idTemp = cell.Edge(position)->Id();
                GenericEdge &edge = *meshPointer->Edge(idTemp);
                if (edge.Point(0)->Id() == id_first_opposite_point && edge.Point(1)->Id() == id_second_opposite_point)
                {
                    third_newFace.AddEdge(&edge);
                    edge.AddFace(&third_newFace);
                    break;
                }
                else if (edge.Point(1)->Id() == id_first_opposite_point && edge.Point(0)->Id() == id_second_opposite_point)
                {
                    third_newFace.AddEdge(&edge);
                    edge.AddFace(&third_newFace);
                    break;
                }
            }

            //****************************************************************************************************
            //Da qui si creano i due tetraedri figli

            cell.SetState(false);
            cell.InitializeChilds(2);

            //Si crea il primo figlio e si settano alcuni suoi membri
            GenericCell &newCell_1 = *meshPointer->CreateCell();
            meshPointer->AddCell(&newCell_1);
            cell.AddChild(&newCell_1);
            newCell_1.SetFather(&cell);
            newCell_1.InheritPropertiesByFather();
            newCell_1.InitializeFaces(4);
            newCell_1.InitializeCells(4);
            newCell_1.InitializePoints(4);
            newCell_1.InitializeEdges(6);


            //Si crea il secondo figlio e si settano alcuni suoi membri
            GenericCell &newCell_2 = *meshPointer->CreateCell();
            meshPointer->AddCell(&newCell_2);
            cell.AddChild(&newCell_2);
            newCell_2.SetFather(&cell);
            newCell_2.InheritPropertiesByFather();
            newCell_2.InitializeFaces(4);
            newCell_2.InitializeCells(4);
            newCell_2.InitializePoints(4);
            newCell_2.InitializeEdges(6);


            //Si aggiungono la facce e viceversa
            newCell_1.AddFace(&third_newFace);
            newCell_2.AddFace(&third_newFace);
            third_newFace.AddCell(&newCell_1);
            third_newFace.AddCell(&newCell_2);
            //Si aggiungono le altre due facce coerentemente con il codice precedente (si sa già quali facce appartengono a quali figli:
            // si vuole che 'newCell_1' abbia 'newEdge_1', first_newFace_1', .....)
            newCell_1.AddFace(&first_newFace_1);
            newCell_2.AddFace(&first_newFace_2);
            first_newFace_1.AddCell(&newCell_1);
            first_newFace_2.AddCell(&newCell_2);
            newCell_1.AddFace(&second_newFace_1);
            newCell_2.AddFace(&second_newFace_2);
            second_newFace_1.AddCell(&newCell_1);
            second_newFace_2.AddCell(&newCell_2);
            flag = 0;
            //Si trova la quarta faccia per entrambi (facce mai toccate) e le si inserisce (e viceversa)
            for (unsigned int position = 0; position < cell.NumberOfFaces(); ++position)
            {
                if (flag == 2)
                    break;
                unsigned int idTemp = cell.Face(position)->Id();
                if (idTemp != id_first_face && idTemp != id_second_face)
                {
                    GenericFace &otherFace = *meshPointer->Face(idTemp);
                    for (int vertex = 0; vertex < otherFace.NumberOfPoints(); ++vertex)
                    {
                        unsigned int idTemp_2 = otherFace.Point(vertex)->Id();
                        if (idTemp_2 == id_first_point)
                        {
                            ++flag;
                            newCell_1.AddFace(&otherFace);
                            otherFace.AddCell(&newCell_1);
                            break;
                        }
                        else if (idTemp_2 == id_second_point)
                        {
                            ++flag;
                            newCell_2.AddFace(&otherFace);
                            otherFace.AddCell(&newCell_2);
                            break;
                        }
                    }
                }
            }
            if(flag != 2)
            {
                Output::PrintErrorMessage("The other two faces for the children of the cell with id %d were not found\n", false, cell.Id());
                return Output::GenericError;
            }

            //Si aggiungono i lati ai tetraedri figli (e viceversa) coerentemente con il codice precedente
            newCell_1.AddEdge(&first_newEdge);
            newCell_1.AddEdge(&second_newEdge);
            first_newEdge.AddCell(&newCell_1);
            second_newEdge.AddCell(&newCell_1);
            newCell_2.AddEdge(&first_newEdge);
            newCell_2.AddEdge(&second_newEdge);
            first_newEdge.AddCell(&newCell_2);
            second_newEdge.AddCell(&newCell_2);
            newCell_1.AddEdge(&newEdge_1);
            newCell_2.AddEdge(&newEdge_2);
            newEdge_1.AddCell(&newCell_1);
            newEdge_2.AddCell(&newCell_2);
            //A partire da codice precedente
            GenericEdge* otherEdge = meshPointer->Edge(first_newFace_1.Edge(1)->Id());
            newCell_1.AddEdge(otherEdge);
            otherEdge->AddCell(&newCell_1);
            otherEdge = meshPointer->Edge(first_newFace_2.Edge(1)->Id());
            newCell_2.AddEdge(otherEdge);
            otherEdge->AddCell(&newCell_2);
            otherEdge = meshPointer->Edge(second_newFace_1.Edge(1)->Id());
            newCell_1.AddEdge(otherEdge);
            otherEdge->AddCell(&newCell_1);
            otherEdge = meshPointer->Edge(second_newFace_2.Edge(1)->Id());
            newCell_2.AddEdge(otherEdge);
            otherEdge->AddCell(&newCell_2);
            otherEdge = meshPointer->Edge(third_newFace.Edge(2)->Id());
            newCell_1.AddEdge(otherEdge);
            otherEdge->AddCell(&newCell_1);
            newCell_2.AddEdge(otherEdge);
            otherEdge->AddCell(&newCell_2);

            //Si aggiungono i punti ai tetraedri e viceversa
            newCell_1.AddPoint(&newPoint);
            newCell_2.AddPoint(&newPoint);
            newPoint.AddCell(&newCell_1);
            newPoint.AddCell(&newCell_2);
            newCell_1.AddPoint(&first_point);
            first_point.AddCell(&newCell_1);
            newCell_2.AddPoint(&second_point);
            second_point.AddCell(&newCell_2);
            newCell_1.AddPoint(&first_opposite_point);
            newCell_2.AddPoint(&first_opposite_point);
            first_opposite_point.AddCell(&newCell_1);
            first_opposite_point.AddCell(&newCell_2);
            newCell_1.AddPoint(&second_opposite_point);
            newCell_2.AddPoint(&second_opposite_point);
            second_opposite_point.AddCell(&newCell_1);
            second_opposite_point.AddCell(&newCell_2);

            //Si aggiungono i tetraedri vicini (ossia quelli con in comune una faccia, non un lato) del padre ai figli coerentemente con la costruzione dei figli
            //Per i tetraedri vicini che hanno il lato più lungo in comune, non si aggiungono qui (questo caso di lato più lungo attivo è contemplato in recupero conformità),
            //ma si aggiungono in recupero conformità
            //Per gli altri tetraedri, se conformi (ossia non tagliati lungo la faccia che in comune qui), si aggiunge figlio che contiene la faccia e
            //i figli del tetraedro che si taglia si aggiungono a quest'ultimo al posto della cella che si sta tagliando;
            //se non conformi, non si aggiungono (si aggiungeranno in conformità)
            GenericFace* face = meshPointer->Face(newCell_1.Face(3)->Id());
            //Con questo 'if' sai se il vicino è stato tagliato lungo questa faccia
            flag = 0;
            if (face->IsActive())
            {
                unsigned int flag_2 = 0;
                for (int position = 0; position < face->NumberOfCells() && flag_2 == 0; ++position)
                {
                    if (face->Cell(position) != NULL)
                    {
                        unsigned int idTemp = face->Cell(position)->Id();
                        if (idTemp != newCell_1.Id() && idTemp != cell.Id())
                        {
                            GenericCell* cell_neig = meshPointer->Cell(idTemp);
                            while (cell_neig->HasChilds())
                            {
                                unsigned int idChild = cell_neig->Child(0)->Id();
                                GenericCell &child = *meshPointer->Cell(idChild);
                                for (int faces = 0; faces < child.NumberOfFaces(); ++faces)
                                {
                                    if (child.Face(faces)->Id() == face->Id())
                                    {
                                        ++flag;
                                        break;
                                    }
                                }
                                if (flag == 0)
                                {
                                    idChild = cell_neig->Child(1)->Id();
                                    cell_neig = meshPointer->Cell(idChild);
                                }
                                else
                                {
                                    cell_neig = meshPointer->Cell(idChild);
                                }
                                flag = 0;
                            }
                            /*flag = 0;
                            //Controllo per controllare se è già dentro a vicini (può succedere: se faccia è in comune con padre e figli e
                            //padre e figlio già inserito in iterazioni precedenti)
                            for (int faces = 0; faces < newCell_1.NumberOfFaces(); ++faces)
                            {
                                unsigned int idTemp_2 = newCell_1.Face(faces)->Id();
                                if (idTemp_2 == cell_neig.Id())
                                {
                                    ++flag;
                                    break;
                                }
                            }*/
                            //if (flag == 0)
                            //{
                            ++flag_2;
                            newCell_1.AddCell(cell_neig);
                            for (int cel = 0; cel < cell_neig->NumberOfCells(); cel++)
                            {
                                if (cell_neig->Cell(cel) != NULL)
                                {
                                    idTemp = cell_neig->Cell(cel)->Id();
                                    if (idTemp == cell.Id())
                                    {
                                        cell_neig->InsertCell(&newCell_1, cel);
                                    }
                                }
                            }
                            //}
                            flag = 0;
                        }
                    }
                }
            }
            face = meshPointer->Face(newCell_2.Face(3)->Id());
            if (face->IsActive())
            {
                unsigned int flag_2 = 0;
                for (int position = 0; position < face->NumberOfCells() && flag_2 == 0; ++position)
                {
                    if (face->Cell(position) != NULL)
                    {
                        unsigned int idTemp = face->Cell(position)->Id();
                        if (idTemp != newCell_2.Id() && idTemp != cell.Id())
                        {
                            GenericCell* cell_neig = meshPointer->Cell(idTemp);
                            while (cell_neig->HasChilds())
                            {
                                unsigned int idChild = cell_neig->Child(0)->Id();
                                GenericCell &child = *meshPointer->Cell(idChild);
                                for (int faces = 0; faces < child.NumberOfFaces(); ++faces)
                                {
                                    if (child.Face(faces)->Id() == face->Id())
                                    {
                                        ++flag;
                                        break;
                                    }
                                }
                                if (flag == 0)
                                {
                                    idChild = cell_neig->Child(1)->Id();
                                    cell_neig = meshPointer->Cell(idChild);
                                }
                                else
                                {
                                    cell_neig = meshPointer->Cell(idChild);
                                }
                                flag = 0;
                            }
                            /*flag = 0;
                            //Controllo per controllare se è già dentro a vicini (può succedere)
                            for (int faces = 0; faces < newCell_2.NumberOfFaces(); ++faces)
                            {
                                unsigned int idTemp_2 = newCell_2.Face(faces)->Id();
                                if (idTemp_2 == cell_neig.Id())
                                {
                                    ++flag;
                                    break;
                                }
                            }
                            if (flag == 0)
                            {*/
                            ++flag_2;
                            newCell_2.AddCell(cell_neig);
                            for (int cel = 0; cel < cell_neig->NumberOfCells(); cel++)
                            {
                                if (cell_neig->Cell(cel) != NULL)
                                {
                                    idTemp = cell_neig->Cell(cel)->Id();
                                    if (idTemp == cell.Id())
                                    {
                                        cell_neig->InsertCell(&newCell_2, cel);
                                    }
                                }
                            }
                            //}
                            flag = 0;
                        }
                    }
                }
            }
            newCell_1.AddCell(&newCell_2);
            newCell_2.AddCell(&newCell_1);

        }
        else
            //****************************************************************************************************
            //caso di lato più lungo tagliato
        {
            to_add = false;

            GenericEdge &newEdge_1 = *meshPointer->Edge(LongestEdge.Child(0)->Id());
            GenericEdge &newEdge_2 = *meshPointer->Edge(LongestEdge.Child(1)->Id());
            GenericPoint &newPoint = *meshPointer->Point(newEdge_1.Point(0)->Id());

            if (first_face.IsActive() && second_face.IsActive())
            {
                //****************************************************************************************************
                //Si taglia la prima faccia
                first_face.SetState(false);
                first_face.InitializeChilds(2);

                //Si creano facce figlie della faccia che si sta tagliando
                GenericFace &first_newFace_1 = *(meshPointer->CreateFace());
                meshPointer->AddFace(&first_newFace_1);
                first_face.AddChild(&first_newFace_1);
                first_newFace_1.SetFather(&first_face);
                first_newFace_1.InheritPropertiesByFather();
                first_newFace_1.AllocateCells(2);
                GenericFace &first_newFace_2 = *(meshPointer->CreateFace());
                meshPointer->AddFace(&first_newFace_2);
                first_face.AddChild(&first_newFace_2);
                first_newFace_2.SetFather(&first_face);
                first_newFace_2.InheritPropertiesByFather();
                first_newFace_2.AllocateCells(2);

                //si aggiungono i lati alle facce e le nuove facce ai lati
                first_newFace_1.InitializeEdges(first_face.NumberOfEdges());
                first_newFace_2.InitializeEdges(first_face.NumberOfEdges());
                first_newFace_1.AddEdge(&newEdge_1);
                first_newFace_2.AddEdge(&newEdge_2);
                newEdge_1.AddFace(&first_newFace_1);
                newEdge_2.AddFace(&first_newFace_2);
                for (unsigned int i = 0; i < first_face.NumberOfEdges(); i++)
                {
                    unsigned int idTemp = first_face.Edge(i)->Id();
                    //Controllo su lato per capire se è quello che si sta tagliando (ossia il più lungo)
                    //no controllo su lato se attivo (perché anche se uno può non esserlo, perché tagliato in un tetraedro che condivide solo un lato,
                    //non è di interesse qui)
                    if (idTemp != idLongestEdge)
                    {
                        GenericEdge &edge = *meshPointer->Edge(idTemp);
                        //Controllo su lato per capire a quale faccia deve appartenere (tenendo conto che 'first_newFace_1' avrà 'first_point', ......
                        if (edge.Point(0)->Id() == id_first_point || edge.Point(1)->Id() == id_first_point)
                        {
                            edge.AddFace(&first_newFace_1);
                            first_newFace_1.AddEdge(&edge);
                        }
                        else if (edge.Point(0)->Id() == id_second_point || edge.Point(1)->Id() == id_second_point)
                        {
                            edge.AddFace(&first_newFace_2);
                            first_newFace_2.AddEdge(&edge);
                        }
                        else
                            return Output::GenericError;
                    }
                }

                //Si aggiungono i punti alle facce e viceversa
                first_newFace_1.InitializePoints(3);
                first_newFace_2.InitializePoints(3);
                first_newFace_1.AddPoint(&newPoint);
                first_newFace_1.AddPoint(&first_point);
                //Per come strutturato, so che il terzo vertice (comune ai due figli) di questa faccia ('first_face') è 'first_opposite_point'
                first_newFace_1.AddPoint(&first_opposite_point);
                first_newFace_2.AddPoint(&newPoint);
                first_newFace_2.AddPoint(&second_point);
                first_newFace_2.AddPoint(&first_opposite_point);
                newPoint.AddFace(&first_newFace_1);
                newPoint.AddFace(&first_newFace_2);
                first_point.AddFace(&first_newFace_1);
                second_point.AddFace(&first_newFace_2);
                first_opposite_point.AddFace(&first_newFace_1);
                first_opposite_point.AddFace(&first_newFace_2);


                //Si crea nuovo lato che deve tagliare la faccia 'first_face'
                GenericEdge &first_newEdge = *(meshPointer->CreateEdge());
                meshPointer->AddEdge(&first_newEdge);
                //Aggiungiamo punto medio al nuovo lato e l'altro vertice (so che è 'first_opposite_point') e viceversa
                first_newEdge.AddPoint(&newPoint);
                newPoint.AddEdge(&first_newEdge);
                first_newEdge.AddPoint(&first_opposite_point);
                first_opposite_point.AddEdge(&first_newEdge);

                //check_presence.push_back(false);

                //Si aggiungono le nuove facce al nuovo lato (e viceversa: non si era ancora inserito il nuovo lato perché non ancora creato quando si trattavano
                //le facce figlie)
                first_newEdge.AddFace(&first_newFace_1);
                first_newEdge.AddFace(&first_newFace_2);
                first_newFace_1.AddEdge(&first_newEdge);
                first_newFace_2.AddEdge(&first_newEdge);

                //****************************************************************************************************
                //Da qui, si inizia a tagliare la seconda faccia che condivide il lato più lungo
                //Si sa già qual è: 'second_face'

                second_face.SetState(false);
                second_face.InitializeChilds(2);

                //Si creano facce figlie
                GenericFace &second_newFace_1 = *(meshPointer->CreateFace());
                meshPointer->AddFace(&second_newFace_1);
                second_face.AddChild(&second_newFace_1);
                second_newFace_1.SetFather(&second_face);
                second_newFace_1.InheritPropertiesByFather();
                second_newFace_1.AllocateCells(2);
                //Le celle che condividono una faccia sono solo due ma non si inseriscono ora (il tetraedro figlio di questa cella sarà inserito dopo e il vicino
                //nel recupero della conformità: come per 'first_face')

                GenericFace &second_newFace_2 = *(meshPointer->CreateFace());
                meshPointer->AddFace(&second_newFace_2);
                second_face.AddChild(&second_newFace_2);
                second_newFace_2.SetFather(&second_face);
                second_newFace_2.InheritPropertiesByFather();
                second_newFace_2.AllocateCells(2);

                //si aggiungono i lati alle facce e le nuove facce ai lati
                second_newFace_1.InitializeEdges(3);
                second_newFace_2.InitializeEdges(3);
                second_newFace_1.AddEdge(&newEdge_1);
                second_newFace_2.AddEdge(&newEdge_2);
                newEdge_1.AddFace(&second_newFace_1);
                newEdge_2.AddFace(&second_newFace_2);
                for (unsigned int i = 0; i < second_face.NumberOfEdges(); i++)
                {
                    unsigned int idTemp = second_face.Edge(i)->Id();
                    //Controllo su lato per capire se è quello che si sta tagliando (ossia il più lungo)
                    //no controllo su lato se attivo (perché anche se uno può non esserlo, perché tagliato in un tetraedro che condivide solo un lato,
                    //non è di interesse qui, ma nel recupero della conformità)
                    if (idTemp != idLongestEdge)
                    {
                        GenericEdge &edge = *meshPointer->Edge(idTemp);
                        //Controllo su lato per capire a quale faccia deve appartenere (tenendo conto che 'first_newFace_1' avrà 'first_point', ......
                        if (edge.Point(0)->Id() == id_first_point || edge.Point(1)->Id() == id_first_point)
                        {
                            edge.AddFace(&second_newFace_1);
                            second_newFace_1.AddEdge(&edge);
                        }
                        else if (edge.Point(0)->Id() == id_second_point || edge.Point(1)->Id() == id_second_point)
                        {
                            edge.AddFace(&second_newFace_2);
                            second_newFace_2.AddEdge(&edge);
                        } else
                            return Output::GenericError;
                    }
                }

                //Si aggiungono i punti alle facce e viceversa
                second_newFace_1.InitializePoints(3);
                second_newFace_2.InitializePoints(3);
                second_newFace_1.AddPoint(&newPoint);
                second_newFace_1.AddPoint(&first_point);
                //Per come strutturato, so che il terzo vertice (comune ai due figli) di questa faccia ('second_face') è 'second_opposite_point'
                second_newFace_1.AddPoint(&second_opposite_point);
                second_newFace_2.AddPoint(&newPoint);
                second_newFace_2.AddPoint(&second_point);
                second_newFace_2.AddPoint(&second_opposite_point);
                newPoint.AddFace(&second_newFace_1);
                newPoint.AddFace(&second_newFace_2);
                first_point.AddFace(&second_newFace_1);
                second_point.AddFace(&second_newFace_2);
                second_opposite_point.AddFace(&second_newFace_1);
                second_opposite_point.AddFace(&second_newFace_2);


                //Si crea nuovo lato che deve tagliare la faccia 'second_face'
                GenericEdge &second_newEdge = *(meshPointer->CreateEdge());
                meshPointer->AddEdge(&second_newEdge);
                //Aggiungiamo punto medio al nuovo lato e l'altro vertice (so che è 'first_opposite_point') e viceversa
                second_newEdge.AddPoint(&newPoint);
                newPoint.AddEdge(&second_newEdge);
                second_newEdge.AddPoint(&second_opposite_point);
                second_opposite_point.AddEdge(&second_newEdge);

                //check_presence.push_back(false);

                //Si aggiungono le nuove facce al nuovo lato (e viceversa: non si era ancora inserito il nuovo lato perché non ancora creato quando si trattavano
                //le facce figlie)
                second_newEdge.AddFace(&second_newFace_1);
                second_newEdge.AddFace(&second_newFace_2);
                second_newFace_1.AddEdge(&second_newEdge);
                second_newFace_2.AddEdge(&second_newEdge);

                //****************************************************************************************************
                //Da qui si crea la terza nuova faccia

                GenericFace &third_newFace = *(meshPointer->CreateFace());
                meshPointer->AddFace(&third_newFace);
                third_newFace.AllocateCells(2);
                //Saranno inseriti i tetredri figli

                //Si aggiungono i punti alla faccia (si conoscono già) e viceversa
                third_newFace.InitializePoints(3);
                third_newFace.AddPoint(&newPoint);
                third_newFace.AddPoint(&first_opposite_point);
                third_newFace.AddPoint(&second_opposite_point);
                first_opposite_point.AddFace(&third_newFace);
                second_opposite_point.AddFace(&third_newFace);
                newPoint.AddFace(&third_newFace);

                //Si aggiungono i lati alla terza faccia e viceversa
                third_newFace.InitializeEdges(3);
                third_newFace.AddEdge(&first_newEdge);
                third_newFace.AddEdge(&second_newEdge);
                first_newEdge.AddFace(&third_newFace);
                second_newEdge.AddFace(&third_newFace);
                //Si deve andare a prendere il terzo lato che non è mai stato trovato (quello con 'first_opposite_point' e 'second_opposite_point'
                //come vertici)
                for (int position = 0; position < cell.NumberOfEdges(); ++position)
                {
                    unsigned int idTemp = cell.Edge(position)->Id();
                    GenericEdge &edge = *meshPointer->Edge(idTemp);
                    if (edge.Point(0)->Id() == id_first_opposite_point && edge.Point(1)->Id() == id_second_opposite_point)
                    {
                        third_newFace.AddEdge(&edge);
                        edge.AddFace(&third_newFace);
                        break;
                    }
                    else if (edge.Point(1)->Id() == id_first_opposite_point && edge.Point(0)->Id() == id_second_opposite_point)
                    {
                        third_newFace.AddEdge(&edge);
                        edge.AddFace(&third_newFace);
                        break;
                    }
                }

                //****************************************************************************************************
                //Da qui si creano i due tetraedri figli

                cell.SetState(false);
                cell.InitializeChilds(2);

                //Si crea il primo figlio e si settano alcuni suoi membri
                GenericCell &newCell_1 = *meshPointer->CreateCell();
                meshPointer->AddCell(&newCell_1);
                cell.AddChild(&newCell_1);
                newCell_1.SetFather(&cell);
                newCell_1.InheritPropertiesByFather();
                newCell_1.InitializeFaces(4);
                newCell_1.InitializeCells(4);
                newCell_1.InitializePoints(4);
                newCell_1.InitializeEdges(6);


                //Si crea il secondo figlio e si settano alcuni suoi membri
                GenericCell &newCell_2 = *meshPointer->CreateCell();
                meshPointer->AddCell(&newCell_2);
                cell.AddChild(&newCell_2);
                newCell_2.SetFather(&cell);
                newCell_2.InheritPropertiesByFather();
                newCell_2.InitializeFaces(4);
                newCell_2.InitializeCells(4);
                newCell_2.InitializePoints(4);
                newCell_2.InitializeEdges(6);


                //Si aggiungono la facce e viceversa
                newCell_1.AddFace(&third_newFace);
                newCell_2.AddFace(&third_newFace);
                third_newFace.AddCell(&newCell_1);
                third_newFace.AddCell(&newCell_2);
                //Si aggiungono le altre due facce coerentemente con il codice precedente (si sa già quali facce appartengono a quali figli:
                // si vuole che 'newCell_1' abbia 'newEdge_1', first_newFace_1', .....)
                newCell_1.AddFace(&first_newFace_1);
                newCell_2.AddFace(&first_newFace_2);
                first_newFace_1.AddCell(&newCell_1);
                first_newFace_2.AddCell(&newCell_2);
                newCell_1.AddFace(&second_newFace_1);
                newCell_2.AddFace(&second_newFace_2);
                second_newFace_1.AddCell(&newCell_1);
                second_newFace_2.AddCell(&newCell_2);
                flag = 0;
                //Si trova la quarta faccia per entrambi (facce mai toccate) e le si inserisce
                for (unsigned int position = 0; position < cell.NumberOfFaces(); ++position)
                {
                    if (flag == 2)
                        break;
                    unsigned int idTemp = cell.Face(position)->Id();
                    if (idTemp != id_first_face && idTemp != id_second_face)
                    {
                        GenericFace &otherFace = *meshPointer->Face(idTemp);
                        for (int vertex = 0; vertex < otherFace.NumberOfPoints(); ++vertex)
                        {
                            unsigned int idTemp_2 = otherFace.Point(vertex)->Id();
                            if (idTemp_2 == id_first_point)
                            {
                                ++flag;
                                newCell_1.AddFace(&otherFace);
                                otherFace.AddCell(&newCell_1);
                                break;
                            }
                            else if (idTemp_2 == id_second_point)
                            {
                                ++flag;
                                newCell_2.AddFace(&otherFace);
                                otherFace.AddCell(&newCell_2);
                                break;
                            }
                        }
                    }
                }
                if(flag != 2)
                {
                    Output::PrintErrorMessage("The other two faces for the children of the cell with id %d were not found\n", false, cell.Id());
                    return Output::GenericError;
                }


                //Si aggiungono i lati ai tetraedri figli (e viceversa) coerentemente con il codice precedente
                newCell_1.AddEdge(&first_newEdge);
                newCell_1.AddEdge(&second_newEdge);
                first_newEdge.AddCell(&newCell_1);
                second_newEdge.AddCell(&newCell_1);
                newCell_2.AddEdge(&first_newEdge);
                newCell_2.AddEdge(&second_newEdge);
                first_newEdge.AddCell(&newCell_2);
                second_newEdge.AddCell(&newCell_2);
                newCell_1.AddEdge(&newEdge_1);
                newCell_2.AddEdge(&newEdge_2);
                newEdge_1.AddCell(&newCell_1);
                newEdge_2.AddCell(&newCell_2);
                GenericEdge* otherEdge = meshPointer->Edge(first_newFace_1.Edge(1)->Id());
                newCell_1.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_1);
                otherEdge = meshPointer->Edge(first_newFace_2.Edge(1)->Id());
                newCell_2.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_2);
                otherEdge = meshPointer->Edge(second_newFace_1.Edge(1)->Id());
                newCell_1.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_1);
                otherEdge = meshPointer->Edge(second_newFace_2.Edge(1)->Id());
                newCell_2.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_2);
                otherEdge = meshPointer->Edge(third_newFace.Edge(2)->Id());
                newCell_1.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_1);
                newCell_2.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_2);


                //Si aggiungono i punti ai tetraedri e viceversa (sostituendo il padre nel vettore dei punti che appartenevano già al padre)
                newCell_1.AddPoint(&newPoint);
                newCell_2.AddPoint(&newPoint);
                newPoint.AddCell(&newCell_1);
                newPoint.AddCell(&newCell_2);
                newCell_1.AddPoint(&first_point);
                first_point.AddCell(&newCell_1);
                newCell_2.AddPoint(&second_point);
                second_point.AddCell(&newCell_2);
                newCell_1.AddPoint(&first_opposite_point);
                newCell_2.AddPoint(&first_opposite_point);
                first_opposite_point.AddCell(&newCell_2);
                first_opposite_point.AddCell(&newCell_1);
                newCell_1.AddPoint(&second_opposite_point);
                newCell_2.AddPoint(&second_opposite_point);
                second_opposite_point.AddCell(&newCell_1);
                second_opposite_point.AddCell(&newCell_2);

                //Si aggiungono i tetraedri vicini (ossia quelli con in comune una faccia, non un lato) del padre ai figli coerentemente con la costruzione dei figli
                //Per i tetraedri vicini che hanno il lato più lungo in comune, non si aggiungono qui (questo caso di lato più lungo attivo è contemplato in recupero conformità),
                //ma si aggiungono in recupero conformità
                //Per gli altri tetraedri, se conformi (ossia non tagliati lungo la faccia che in comune qui), si aggiunge figlio che contiene la faccia e
                //i figli del tetraedro che si taglia si aggiungono a quest'ultimo al posto della cella che si sta taagliando;
                //se non conformi, non si aggiungono (si aggiungeranno in conformità)
                GenericFace* face = meshPointer->Face(newCell_1.Face(3)->Id());
                //Con questo 'if' sai se il vicino è stato tagliato lungo questa faccia
                flag = 0;
                if (face->IsActive())
                {
                    unsigned int flag_2 = 0;
                    for (int position = 0; position < face->NumberOfCells() && flag_2 == 0; ++position)
                    {
                        if(face->Cell(position) != NULL)
                        {
                            unsigned int idTemp = face->Cell(position)->Id();
                            if (idTemp != newCell_1.Id() && idTemp != cell.Id())
                            {
                                GenericCell* cell_neig = meshPointer->Cell(idTemp);
                                while (cell_neig->HasChilds())
                                {
                                    unsigned int idChild = cell_neig->Child(0)->Id();
                                    GenericCell &child = *meshPointer->Cell(idChild);
                                    for (int faces = 0; faces < child.NumberOfFaces(); ++faces)
                                    {
                                        if (child.Face(faces)->Id() == face->Id())
                                        {
                                            ++flag;
                                            break;
                                        }
                                    }
                                    if (flag == 0)
                                    {
                                        idChild = cell_neig->Child(1)->Id();
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    else
                                    {
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    flag = 0;
                                }
                                /*flag = 0;
                                //Controllo per controllare se è già dentro a vicini (può succedere)
                                for (int faces = 0; faces < newCell_1.NumberOfFaces(); ++faces)
                                {
                                    unsigned int idTemp_2 = newCell_1.Face(faces)->Id();
                                    if (idTemp_2 == cell_neig.Id())
                                    {
                                        ++flag;
                                        break;
                                    }
                                }
                                if (flag == 0)
                                {*/
                                ++flag_2;
                                newCell_1.AddCell(cell_neig);
                                for (int cel = 0; cel < cell_neig->NumberOfCells(); cel++)
                                {
                                    if(cell_neig->Cell(cel) != NULL)
                                    {
                                        idTemp = cell_neig->Cell(cel)->Id();
                                        if (idTemp == cell.Id())
                                        {
                                            cell_neig->InsertCell(&newCell_1, cel);
                                        }
                                    }
                                }
                                //}
                                flag = 0;
                            }
                        }
                    }
                }
                face = meshPointer->Face(newCell_2.Face(3)->Id());
                if (face->IsActive())
                {
                    unsigned int flag_2 = 0;
                    for (int position = 0; position < face->NumberOfCells() && flag_2 == 0; ++position)
                    {
                        if(face->Cell(position) != NULL)
                        {
                            unsigned int idTemp = face->Cell(position)->Id();
                            if (idTemp != newCell_2.Id() && idTemp != cell.Id())
                            {
                                GenericCell* cell_neig = meshPointer->Cell(idTemp);
                                while (cell_neig->HasChilds())
                                {
                                    unsigned int idChild = cell_neig->Child(0)->Id();
                                    GenericCell &child = *meshPointer->Cell(idChild);
                                    for (int faces = 0; faces < child.NumberOfFaces(); ++faces)
                                    {
                                        if (child.Face(faces)->Id() == face->Id())
                                        {
                                            ++flag;
                                            break;
                                        }
                                    }
                                    if (flag == 0)
                                    {
                                        idChild = cell_neig->Child(1)->Id();
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    else
                                    {
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    flag = 0;
                                }
                                /*flag = 0;
                                //Controllo per controllare se è già dentro a vicini (può succedere)
                                for (int faces = 0; faces < newCell_2.NumberOfFaces(); ++faces)
                                {
                                    unsigned int idTemp_2 = newCell_2.Face(faces)->Id();
                                    if (idTemp_2 == cell_neig.Id())
                                    {
                                        ++flag;
                                        break;
                                    }
                                }
                                if (flag == 0) {*/
                                ++flag_2;
                                newCell_2.AddCell(cell_neig);
                                for (int cel = 0; cel < cell_neig->NumberOfCells(); cel++)
                                {
                                    if(cell_neig->Cell(cel) != NULL)
                                    {
                                        idTemp = cell_neig->Cell(cel)->Id();
                                        if (idTemp == cell.Id())
                                        {
                                            cell_neig->InsertCell(&newCell_2, cel);
                                        }
                                    }
                                }
                                //}
                                flag = 0;
                            }
                        }
                    }
                }
                newCell_1.AddCell(&newCell_2);
                newCell_2.AddCell(&newCell_1);
            }
                //****************************************************************************************************
                //Altro caso
            else if (first_face.IsActive() && !second_face.IsActive())
            {
                //Si taglia la prima faccia
                first_face.SetState(false);
                first_face.InitializeChilds(2);

                //Si creano facce figlie della faccia che si sta tagliando
                GenericFace &first_newFace_1 = *(meshPointer->CreateFace());
                meshPointer->AddFace(&first_newFace_1);
                first_face.AddChild(&first_newFace_1);
                first_newFace_1.SetFather(&first_face);
                first_newFace_1.InheritPropertiesByFather();
                first_newFace_1.AllocateCells(2);
                GenericFace &first_newFace_2 = *(meshPointer->CreateFace());
                meshPointer->AddFace(&first_newFace_2);
                first_face.AddChild(&first_newFace_2);
                first_newFace_2.SetFather(&first_face);
                first_newFace_2.InheritPropertiesByFather();
                first_newFace_2.AllocateCells(2);

                //si aggiungono i lati alle facce e le nuove facce ai lati
                first_newFace_1.InitializeEdges(first_face.NumberOfEdges());
                first_newFace_2.InitializeEdges(first_face.NumberOfEdges());
                first_newFace_1.AddEdge(&newEdge_1);
                first_newFace_2.AddEdge(&newEdge_2);
                for (unsigned int i = 0; i < first_face.NumberOfEdges(); i++)
                {
                    unsigned int idTemp = first_face.Edge(i)->Id();
                    //Controllo su lato per capire se è quello che si sta tagliando (ossia il più lungo)
                    //no controllo su lato se attivo (perché anche se uno può non esserlo, perché tagliato in un tetraedro che condivide solo un lato,
                    //non è di interesse qui)
                    if (idTemp != idLongestEdge)
                    {
                        GenericEdge &edge = *meshPointer->Edge(idTemp);
                        //Controllo su lato per capire a quale faccia deve appartenere (tenendo conto che 'first_newFace_1' avrà 'first_point', ......
                        if (edge.Point(0)->Id() == id_first_point || edge.Point(1)->Id() == id_first_point)
                        {
                            edge.AddFace(&first_newFace_1);
                            first_newFace_1.AddEdge(&edge);
                        }
                        else if (edge.Point(0)->Id() == id_second_point || edge.Point(1)->Id() == id_second_point)
                        {
                            edge.AddFace(&first_newFace_2);
                            first_newFace_2.AddEdge(&edge);
                        }
                        else
                            return Output::GenericError;
                    }
                }
                //Prima non si era inserito la faccia massima perché sarebbe stata tagliata
                newEdge_1.AddFace(&first_newFace_1);
                newEdge_2.AddFace(&first_newFace_2);

                //Si aggiungono i punti alle facce e viceversa
                first_newFace_1.InitializePoints(3);
                first_newFace_2.InitializePoints(3);
                first_newFace_1.AddPoint(&newPoint);
                first_newFace_1.AddPoint(&first_point);
                //Per come strutturato, so che il terzo vertice (comune ai due figli) di questa faccia ('first_face') è 'first_opposite_point'
                first_newFace_1.AddPoint(&first_opposite_point);
                first_newFace_2.AddPoint(&newPoint);
                first_newFace_2.AddPoint(&second_point);
                first_newFace_2.AddPoint(&first_opposite_point);
                newPoint.AddFace(&first_newFace_1);
                newPoint.AddFace(&first_newFace_2);
                first_point.AddFace(&first_newFace_1);
                second_point.AddFace(&first_newFace_2);
                first_opposite_point.AddFace(&first_newFace_1);
                first_opposite_point.AddFace(&first_newFace_2);


                //Si crea nuovo lato che deve tagliare la faccia 'first_face'
                GenericEdge &first_newEdge = *(meshPointer->CreateEdge());
                meshPointer->AddEdge(&first_newEdge);
                //Aggiungiamo punto medio al nuovo lato e l'altro vertice (so che è 'first_opposite_point') e viceversa
                first_newEdge.AddPoint(&newPoint);
                newPoint.AddEdge(&first_newEdge);
                first_newEdge.AddPoint(&first_opposite_point);
                first_opposite_point.AddEdge(&first_newEdge);

                //check_presence.push_back(false);

                //Si aggiungono le nuove facce al nuovo lato (e viceversa: non si era ancora inserito il nuovo lato perché non ancora creato quando si trattavano
                //le facce figlie)
                first_newEdge.AddFace(&first_newFace_1);
                first_newEdge.AddFace(&first_newFace_2);
                first_newFace_1.AddEdge(&first_newEdge);
                first_newFace_2.AddEdge(&first_newEdge);


                //****************************************************************************************************
                //Da qui, si prende la seconda faccia che condivide il lato più lungo (già tagliata)
                //Si sa già qual è: 'second_face'

                unsigned int id_second_nF = second_face.Child(0)->Id();
                GenericFace &second_newFace_1 = *(meshPointer->Face(id_second_nF));
                id_second_nF = second_face.Child(1)->Id();
                GenericFace &second_newFace_2 = *(meshPointer->Face(id_second_nF));

                //Si prende nuovo lato che taglia 'second_face'
                id_second_nF = second_newFace_1.Edge(2)->Id();
                GenericEdge &second_newEdge = *(meshPointer->Edge(id_second_nF));

                //****************************************************************************************************
                //Da qui si crea la terza nuova faccia

                GenericFace &third_newFace = *(meshPointer->CreateFace());
                meshPointer->AddFace(&third_newFace);
                third_newFace.AllocateCells(2);
                //Saranno inseriti i tetredri figli

                //Si aggiungono i punti alla faccia (si conoscono già) e viceversa
                third_newFace.InitializePoints(3);
                third_newFace.AddPoint(&newPoint);
                third_newFace.AddPoint(&first_opposite_point);
                third_newFace.AddPoint(&second_opposite_point);
                first_opposite_point.AddFace(&third_newFace);
                second_opposite_point.AddFace(&third_newFace);
                newPoint.AddFace(&third_newFace);

                //Si aggiungono i lati alla terza faccia e viceversa
                third_newFace.InitializeEdges(3);
                third_newFace.AddEdge(&first_newEdge);
                third_newFace.AddEdge(&second_newEdge);
                first_newEdge.AddFace(&third_newFace);
                second_newEdge.AddFace(&third_newFace);
                //Si deve andare a prendere il terzo lato che non è mai stato trovato (quello con 'first_opposite_point' e 'second_opposite_point'
                //come vertici)
                for (int position = 0; position < cell.NumberOfEdges(); ++position)
                {
                    unsigned int idTemp = cell.Edge(position)->Id();
                    GenericEdge &edge = *meshPointer->Edge(idTemp);
                    if (edge.Point(0)->Id() == id_first_opposite_point && edge.Point(1)->Id() == id_second_opposite_point)
                    {
                        third_newFace.AddEdge(&edge);
                        edge.AddFace(&third_newFace);
                        break;
                    }
                    else if (edge.Point(1)->Id() == id_first_opposite_point && edge.Point(0)->Id() == id_second_opposite_point)
                    {
                        third_newFace.AddEdge(&edge);
                        edge.AddFace(&third_newFace);
                        break;
                    }
                }

                //****************************************************************************************************
                //Da qui si creano i due tetraedri figli

                cell.SetState(false);
                cell.InitializeChilds(2);

                //Si crea il primo figlio e si settano alcuni suoi membri
                GenericCell &newCell_1 = *meshPointer->CreateCell();
                meshPointer->AddCell(&newCell_1);
                cell.AddChild(&newCell_1);
                newCell_1.SetFather(&cell);
                newCell_1.InheritPropertiesByFather();
                newCell_1.InitializeFaces(4);
                newCell_1.InitializeCells(4);
                newCell_1.InitializePoints(4);
                newCell_1.InitializeEdges(6);


                //Si crea il secondo figlio e si settano alcuni suoi membri
                GenericCell &newCell_2 = *meshPointer->CreateCell();
                meshPointer->AddCell(&newCell_2);
                cell.AddChild(&newCell_2);
                newCell_2.SetFather(&cell);
                newCell_2.InheritPropertiesByFather();
                newCell_2.InitializeFaces(4);
                newCell_2.InitializeCells(4);
                newCell_2.InitializePoints(4);
                newCell_2.InitializeEdges(6);


                //Si aggiungono la facce e viceversa
                newCell_1.AddFace(&third_newFace);
                newCell_2.AddFace(&third_newFace);
                third_newFace.AddCell(&newCell_1);
                third_newFace.AddCell(&newCell_2);
                //Si aggiungono le altre due facce coerentemente con il codice precedente (si sa già quali facce appartengono a quali figli:
                //si vuole che 'newCell_1' abbia 'newEdge_1', first_newFace_1', .....)
                newCell_1.AddFace(&first_newFace_1);
                newCell_2.AddFace(&first_newFace_2);
                first_newFace_1.AddCell(&newCell_1);
                first_newFace_2.AddCell(&newCell_2);
                newCell_1.AddFace(&second_newFace_1);
                newCell_2.AddFace(&second_newFace_2);
                second_newFace_1.AddCell(&newCell_1);
                second_newFace_2.AddCell(&newCell_2);
                flag = 0;
                //Si trova la quarta faccia per entrambi (facce mai toccate) e le si inserisce (e viceversa)
                for (unsigned int position = 0; position < cell.NumberOfFaces(); ++position)
                {
                    if (flag == 2)
                        break;
                    unsigned int idTemp = cell.Face(position)->Id();
                    if (idTemp != id_first_face && idTemp != id_second_face)
                    {
                        GenericFace &otherFace = *meshPointer->Face(idTemp);
                        for (int vertex = 0; vertex < otherFace.NumberOfPoints(); ++vertex)
                        {
                            unsigned int idTemp_2 = otherFace.Point(vertex)->Id();
                            if (idTemp_2 == id_first_point)
                            {
                                ++flag;
                                newCell_1.AddFace(&otherFace);
                                otherFace.AddCell(&newCell_1);
                                break;
                            }
                            else if (idTemp_2 == id_second_point)
                            {
                                ++flag;
                                newCell_2.AddFace(&otherFace);
                                otherFace.AddCell(&newCell_2);
                                break;
                            }
                        }
                    }
                }
                if(flag != 2)
                {
                    Output::PrintErrorMessage("The other two faces for the children of the cell with id %d were not found\n", false, cell.Id());
                    return Output::GenericError;
                }


                //Si aggiungono i lati ai tetraedri figli (e viceversa) coerentemente con il codice precedente
                newCell_1.AddEdge(&first_newEdge);
                newCell_1.AddEdge(&second_newEdge);
                first_newEdge.AddCell(&newCell_1);
                second_newEdge.AddCell(&newCell_1);
                newCell_2.AddEdge(&first_newEdge);
                newCell_2.AddEdge(&second_newEdge);
                first_newEdge.AddCell(&newCell_2);
                second_newEdge.AddCell(&newCell_2);
                newCell_1.AddEdge(&newEdge_1);
                newCell_2.AddEdge(&newEdge_2);
                newEdge_1.AddCell(&newCell_1);
                newEdge_2.AddCell(&newCell_2);
                //A partire da codice precedente
                GenericEdge* otherEdge = meshPointer->Edge(first_newFace_1.Edge(1)->Id());
                newCell_1.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_1);
                otherEdge = meshPointer->Edge(first_newFace_2.Edge(1)->Id());
                newCell_2.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_2);
                otherEdge = meshPointer->Edge(second_newFace_1.Edge(1)->Id());
                newCell_1.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_1);
                otherEdge = meshPointer->Edge(second_newFace_2.Edge(1)->Id());
                newCell_2.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_2);
                otherEdge = meshPointer->Edge(third_newFace.Edge(2)->Id());
                newCell_1.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_1);
                newCell_2.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_2);


                //Si aggiungono i punti ai tetraedri e viceversa
                newCell_1.AddPoint(&newPoint);
                newCell_2.AddPoint(&newPoint);
                newPoint.AddCell(&newCell_1);
                newPoint.AddCell(&newCell_2);
                newCell_1.AddPoint(&first_point);
                first_point.AddCell(&newCell_1);
                newCell_2.AddPoint(&second_point);
                second_point.AddCell(&newCell_2);
                newCell_1.AddPoint(&first_opposite_point);
                newCell_2.AddPoint(&first_opposite_point);
                first_opposite_point.AddCell(&newCell_2);
                first_opposite_point.AddCell(&newCell_1);
                newCell_1.AddPoint(&second_opposite_point);
                newCell_2.AddPoint(&second_opposite_point);
                second_opposite_point.AddCell(&newCell_1);
                second_opposite_point.AddCell(&newCell_2);

                //Si aggiungono i tetraedri vicini del padre di (ancora conformi ai figli) ai figli coerentemente con la costruzione dei figli
                //Se la faccia (non con il lato più lungo della cella padre) che è in comune con il tetraedro vicino è stata tagliata nel tetraedro vicino,
                //allora non si mette, si metterà in recupero conformità
                //Se non è stata tagliata, allora tra i vicini si mette il tetraedo (figlio o no) che contiene quella faccia
                //Per quanto riguarda 'second_face', si mettono i due vicini conformi
                GenericFace* face = meshPointer->Face(newCell_1.Face(3)->Id());
                //Con questo 'if' sai se il vicino è stato tagliato lungo questa faccia
                flag = 0;
                if (face->IsActive())
                {
                    unsigned int flag_2 = 0;
                    for (int position = 0; position < face->NumberOfCells() && flag_2 == 0; ++position)
                    {
                        if(face->Cell(position) != NULL)
                        {
                            unsigned int idTemp = face->Cell(position)->Id();
                            if (idTemp != newCell_1.Id() && idTemp != cell.Id())
                            {
                                GenericCell* cell_neig = meshPointer->Cell(idTemp);
                                while (cell_neig->HasChilds())
                                {
                                    unsigned int idChild = cell_neig->Child(0)->Id();
                                    GenericCell &child = *meshPointer->Cell(idChild);
                                    for (int faces = 0; faces < child.NumberOfFaces(); ++faces)
                                    {
                                        if (child.Face(faces)->Id() == face->Id())
                                        {
                                            ++flag;
                                            break;
                                        }
                                    }
                                    if (flag == 0)
                                    {
                                        idChild = cell_neig->Child(1)->Id();
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    else
                                    {
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    flag = 0;
                                }
                                /*flag = 0;
                                //Controllo per controllare se è già dentro a vicini (può succedere)
                                for (int faces = 0; faces < newCell_1.NumberOfFaces(); ++faces)
                                {
                                    unsigned int idTemp_2 = newCell_1.Face(faces)->Id();
                                    if (idTemp_2 == cell_neig.Id())
                                    {
                                        ++flag;
                                        break;
                                    }
                                }
                                if (flag == 0)
                                {*/
                                ++flag_2;
                                newCell_1.AddCell(cell_neig);
                                for (int cel = 0; cel < cell_neig->NumberOfCells(); cel++)
                                {
                                    if(cell_neig->Cell(cel) != NULL)
                                    {
                                        idTemp = cell_neig->Cell(cel)->Id();
                                        if (idTemp == cell.Id())
                                        {
                                            cell_neig->InsertCell(&newCell_1, cel);
                                        }
                                    }
                                }
                                //}
                                flag = 0;
                            }
                        }
                    }
                }
                face = meshPointer->Face(newCell_2.Face(3)->Id());
                if (face->IsActive())
                {
                    unsigned int flag_2 = 0;
                    for (int position = 0; position < face->NumberOfCells() && flag_2 == 0; ++position)
                    {
                        if(face->Cell(position) != NULL)
                        {
                            unsigned int idTemp = face->Cell(position)->Id();
                            if (idTemp != newCell_2.Id() && idTemp != cell.Id()) {
                                GenericCell* cell_neig = meshPointer->Cell(idTemp);
                                while (cell_neig->HasChilds())
                                {
                                    unsigned int idChild = cell_neig->Child(0)->Id();
                                    GenericCell &child = *meshPointer->Cell(idChild);
                                    for (int faces = 0; faces < child.NumberOfFaces(); ++faces)
                                    {
                                        if (child.Face(faces)->Id() == face->Id())
                                        {
                                            ++flag;
                                            break;
                                        }
                                    }
                                    if (flag == 0)
                                    {
                                        idChild = cell_neig->Child(1)->Id();
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    else
                                    {
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    flag = 0;
                                }
                                /*flag = 0;
                                //Controllo per controllare se è già dentro a vicini (può succedere)
                                for (int faces = 0; faces < newCell_2.NumberOfFaces(); ++faces)
                                {
                                    unsigned int idTemp_2 = newCell_2.Face(faces)->Id();
                                    if (idTemp_2 == cell_neig.Id())
                                    {
                                        ++flag;
                                        break;
                                    }
                                }
                                if (flag == 0)
                                {*/
                                ++flag_2;
                                newCell_2.AddCell(cell_neig);
                                for (int cel = 0; cel < cell_neig->NumberOfCells(); cel++)
                                {
                                    if(cell_neig->Cell(cel) != NULL)
                                    {
                                        idTemp = cell_neig->Cell(cel)->Id();
                                        if (idTemp == cell.Id())
                                        {
                                            cell_neig->InsertCell(&newCell_2, cel);
                                        }
                                    }
                                }
                                //}
                                flag = 0;
                            }
                        }
                    }
                }
                //Si è recuperata la conformità
                GenericCell* cell_neig = meshPointer->Cell(first_newEdge.Cell(0)->Id());
                newCell_1.AddCell(cell_neig);
                cell_neig->AddCell(&newCell_1);
                cell_neig = meshPointer->Cell(first_newEdge.Cell(1)->Id());
                newCell_2.AddCell(cell_neig);
                cell_neig->AddCell(&newCell_2);
                newCell_1.AddCell(&newCell_2);
                newCell_2.AddCell(&newCell_1);

            }
            else if (!first_face.IsActive() && second_face.IsActive())
                //****************************************************************************************************
                //Altro caso
            {
                //Si prendono i figli della prima faccia
                unsigned int id_first_nF = first_face.Child(0)->Id();
                GenericFace &first_newFace_1 = *(meshPointer->Face(id_first_nF));
                id_first_nF = first_face.Child(1)->Id();
                GenericFace &first_newFace_2 = *(meshPointer->Face(id_first_nF));

                //Si prende nuovo lato che taglia 'first_face'
                id_first_nF = first_newFace_1.Edge(2)->Id();
                GenericEdge &first_newEdge = *(meshPointer->Edge(id_first_nF));

                //****************************************************************************************************
                //Da qui, si inizia a tagliare la seconda faccia che condivide il lato più lungo
                //Si sa già qual è: 'second_face'

                second_face.SetState(false);
                second_face.InitializeChilds(2);

                //Si creano facce figlie
                GenericFace &second_newFace_1 = *(meshPointer->CreateFace());
                meshPointer->AddFace(&second_newFace_1);
                second_face.AddChild(&second_newFace_1);
                second_newFace_1.SetFather(&second_face);
                second_newFace_1.InheritPropertiesByFather();
                second_newFace_1.AllocateCells(2);
                //Le celle che condividono una faccia sono solo due ma non si inseriscono ora (il tetraedro figlio di questa cella sarà inserito dopo e il vicino
                //nel recupero della conformità)

                GenericFace &second_newFace_2 = *(meshPointer->CreateFace());
                meshPointer->AddFace(&second_newFace_2);
                second_face.AddChild(&second_newFace_2);
                second_newFace_2.SetFather(&second_face);
                second_newFace_2.InheritPropertiesByFather();
                second_newFace_2.AllocateCells(2);

                //si aggiungono i lati alle facce e le nuove facce ai lati
                second_newFace_1.InitializeEdges(3);
                second_newFace_2.InitializeEdges(3);
                second_newFace_1.AddEdge(&newEdge_1);
                second_newFace_2.AddEdge(&newEdge_2);
                for (unsigned int i = 0; i < second_face.NumberOfEdges(); i++)
                {
                    unsigned int idTemp = second_face.Edge(i)->Id();
                    //Controllo su lato per capire se è quello che si sta tagliando (ossia il più lungo)
                    //no controllo su lato se attivo (perché anche se uno può non esserlo, perché tagliato in un tetraedro che condivide solo un lato,
                    //non è di interesse qui, ma nel recupero della conformità)
                    if (idTemp != idLongestEdge) {
                        GenericEdge &edge = *meshPointer->Edge(idTemp);
                        //Controllo su lato per capire a quale faccia deve appartenere (tenendo conto che 'first_newFace_1' avrà 'first_point', ......
                        if (edge.Point(0)->Id() == id_first_point || edge.Point(1)->Id() == id_first_point)
                        {
                            edge.AddFace(&second_newFace_1);
                            second_newFace_1.AddEdge(&edge);
                        }
                        else if (edge.Point(0)->Id() == id_second_point || edge.Point(1)->Id() == id_second_point)
                        {
                            edge.AddFace(&second_newFace_2);
                            second_newFace_2.AddEdge(&edge);
                        }
                        else
                            return Output::GenericError;
                    }
                }
                newEdge_1.AddFace(&second_newFace_1);
                newEdge_2.AddFace(&second_newFace_2);

                //Si aggiungono i punti alle facce e viceversa
                second_newFace_1.InitializePoints(3);
                second_newFace_2.InitializePoints(3);
                second_newFace_1.AddPoint(&newPoint);
                second_newFace_1.AddPoint(&first_point);
                //Per come strutturato, so che il terzo vertice (comune ai due figli) di questa faccia ('second_face') è 'second_opposite_point'
                second_newFace_1.AddPoint(&second_opposite_point);
                second_newFace_2.AddPoint(&newPoint);
                second_newFace_2.AddPoint(&second_point);
                second_newFace_2.AddPoint(&second_opposite_point);
                newPoint.AddFace(&second_newFace_1);
                newPoint.AddFace(&second_newFace_2);
                first_point.AddFace(&second_newFace_1);
                second_point.AddFace(&second_newFace_2);
                second_opposite_point.AddFace(&second_newFace_1);
                second_opposite_point.AddFace(&second_newFace_2);


                //Si crea nuovo lato che deve tagliare la faccia 'second_face'
                GenericEdge &second_newEdge = *(meshPointer->CreateEdge());
                meshPointer->AddEdge(&second_newEdge);
                //Aggiungiamo punto medio al nuovo lato e l'altro vertice (so che è 'first_opposite_point') e viceversa
                second_newEdge.AddPoint(&newPoint);
                newPoint.AddEdge(&second_newEdge);
                second_newEdge.AddPoint(&second_opposite_point);
                second_opposite_point.AddEdge(&second_newEdge);

                //check_presence.push_back(false);

                //Si aggiungono le nuove facce al nuovo lato (e viceversa: non si era ancora inserito il nuovo lato perché non ancora creato quando si trattavano
                //le facce figlie)
                second_newEdge.AddFace(&second_newFace_1);
                second_newEdge.AddFace(&second_newFace_2);
                second_newFace_1.AddEdge(&second_newEdge);
                second_newFace_2.AddEdge(&second_newEdge);

                //****************************************************************************************************
                //Da qui si crea la terza nuova faccia

                GenericFace &third_newFace = *(meshPointer->CreateFace());
                meshPointer->AddFace(&third_newFace);
                third_newFace.AllocateCells(2);
                //Saranno inseriti i tetredri figli

                //Si aggiungono i punti alla faccia (si conoscono già) e viceversa
                third_newFace.InitializePoints(3);
                third_newFace.AddPoint(&newPoint);
                third_newFace.AddPoint(&first_opposite_point);
                third_newFace.AddPoint(&second_opposite_point);
                first_opposite_point.AddFace(&third_newFace);
                second_opposite_point.AddFace(&third_newFace);
                newPoint.AddFace(&third_newFace);

                //Si aggiungono i lati alla terza faccia e viceversa
                third_newFace.InitializeEdges(3);
                third_newFace.AddEdge(&first_newEdge);
                third_newFace.AddEdge(&second_newEdge);
                first_newEdge.AddFace(&third_newFace);
                second_newEdge.AddFace(&third_newFace);
                //Si deve andare a prendere il terzo lato che non è mai stato trovato (quello con 'first_opposite_point' e 'second_opposite_point'
                //come vertici)
                for (int position = 0; position < cell.NumberOfEdges(); ++position)
                {
                    unsigned int idTemp = cell.Edge(position)->Id();
                    GenericEdge &edge = *meshPointer->Edge(idTemp);
                    if (edge.Point(0)->Id() == id_first_opposite_point &&
                        edge.Point(1)->Id() == id_second_opposite_point)
                    {
                        third_newFace.AddEdge(&edge);
                        edge.AddFace(&third_newFace);
                        break;
                    }
                    else if (edge.Point(1)->Id() == id_first_opposite_point && edge.Point(0)->Id() == id_second_opposite_point)
                    {
                        third_newFace.AddEdge(&edge);
                        edge.AddFace(&third_newFace);
                        break;
                    }
                }


                //****************************************************************************************************
                //Da qui si creano i due tetraedri figli

                cell.SetState(false);
                cell.InitializeChilds(2);

                //Si crea il primo figlio e si settano alcuni suoi membri
                GenericCell &newCell_1 = *meshPointer->CreateCell();
                meshPointer->AddCell(&newCell_1);
                cell.AddChild(&newCell_1);
                newCell_1.SetFather(&cell);
                newCell_1.InheritPropertiesByFather();
                newCell_1.InitializeFaces(4);
                newCell_1.InitializeCells(4);
                newCell_1.InitializePoints(4);
                newCell_1.InitializeEdges(6);


                //Si crea il secondo figlio e si settano alcuni suoi membri
                GenericCell &newCell_2 = *meshPointer->CreateCell();
                meshPointer->AddCell(&newCell_2);
                cell.AddChild(&newCell_2);
                newCell_2.SetFather(&cell);
                newCell_2.InheritPropertiesByFather();
                newCell_2.InitializeFaces(4);
                newCell_2.InitializeCells(4);
                newCell_2.InitializePoints(4);
                newCell_2.InitializeEdges(6);


                //Si aggiungono la facce e viceversa
                newCell_1.AddFace(&third_newFace);
                newCell_2.AddFace(&third_newFace);
                third_newFace.AddCell(&newCell_1);
                third_newFace.AddCell(&newCell_2);
                //Si aggiungono le altre due facce coerentemente con il codice precedente (si sa già quali facce appartengono a quali figli:
                //si vuole che 'newCell_1' abbia 'newEdge_1', first_newFace_1', .....)
                newCell_1.AddFace(&first_newFace_1);
                newCell_2.AddFace(&first_newFace_2);
                first_newFace_1.AddCell(&newCell_1);
                first_newFace_2.AddCell(&newCell_2);
                newCell_1.AddFace(&second_newFace_1);
                newCell_2.AddFace(&second_newFace_2);
                second_newFace_1.AddCell(&newCell_1);
                second_newFace_2.AddCell(&newCell_2);
                flag = 0;
                //Si trova la quarta faccia per entrambi (facce mai toccate) e le si inserisce (e viceversa)
                for (unsigned int position = 0; position < cell.NumberOfFaces(); ++position)
                {
                    if (flag == 2)
                        break;
                    unsigned int idTemp = cell.Face(position)->Id();
                    if (idTemp != id_first_face && idTemp != id_second_face)
                    {
                        GenericFace &otherFace = *meshPointer->Face(idTemp);
                        for (int vertex = 0; vertex < otherFace.NumberOfPoints(); ++vertex)
                        {
                            unsigned int idTemp_2 = otherFace.Point(vertex)->Id();
                            if (idTemp_2 == id_first_point)
                            {
                                ++flag;
                                newCell_1.AddFace(&otherFace);
                                otherFace.AddCell(&newCell_1);
                                break;
                            } else if (idTemp_2 == id_second_point)
                            {
                                ++flag;
                                newCell_2.AddFace(&otherFace);
                                otherFace.AddCell(&newCell_2);
                                break;
                            }
                        }
                    }
                }
                if(flag != 2)
                {
                    Output::PrintErrorMessage("The other two faces for the children of the cell with id %d were not found\n", false, cell.Id());
                    return Output::GenericError;
                }


                //Si aggiungono i lati ai tetraedri figli (e viceversa) coerentemente con il codice precedente
                newCell_1.AddEdge(&first_newEdge);
                newCell_1.AddEdge(&second_newEdge);
                first_newEdge.AddCell(&newCell_1);
                second_newEdge.AddCell(&newCell_1);
                newCell_2.AddEdge(&first_newEdge);
                newCell_2.AddEdge(&second_newEdge);
                first_newEdge.AddCell(&newCell_2);
                second_newEdge.AddCell(&newCell_2);
                newCell_1.AddEdge(&newEdge_1);
                newCell_2.AddEdge(&newEdge_2);
                newEdge_1.AddCell(&newCell_1);
                newEdge_2.AddCell(&newCell_2);
                //A partire da codice precedente
                GenericEdge* otherEdge = meshPointer->Edge(first_newFace_1.Edge(1)->Id());
                newCell_1.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_1);
                otherEdge = meshPointer->Edge(first_newFace_2.Edge(1)->Id());
                newCell_2.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_2);
                otherEdge = meshPointer->Edge(second_newFace_1.Edge(1)->Id());
                newCell_1.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_1);
                otherEdge = meshPointer->Edge(second_newFace_2.Edge(1)->Id());
                newCell_2.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_2);
                otherEdge = meshPointer->Edge(third_newFace.Edge(2)->Id());
                newCell_1.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_1);
                newCell_2.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_2);


                //Si aggiungono i punti ai tetraedri e viceversa
                newCell_1.AddPoint(&newPoint);
                newCell_2.AddPoint(&newPoint);
                newPoint.AddCell(&newCell_1);
                newPoint.AddCell(&newCell_2);
                newCell_1.AddPoint(&first_point);
                first_point.AddCell(&newCell_1);
                newCell_2.AddPoint(&second_point);
                second_point.AddCell(&newCell_2);
                newCell_1.AddPoint(&first_opposite_point);
                newCell_2.AddPoint(&first_opposite_point);
                first_opposite_point.AddCell(&newCell_2);
                first_opposite_point.AddCell(&newCell_1);
                newCell_1.AddPoint(&second_opposite_point);
                newCell_2.AddPoint(&second_opposite_point);
                second_opposite_point.AddCell(&newCell_1);
                second_opposite_point.AddCell(&newCell_2);

                //Si aggiungono i tetraedri vicini del padre (ancora conformi ai figli) ai figli coerentemente con la costruzione dei figli
                //Se la faccia (non con il lato più lungo della cella padre) che è in comune con il tetraedro vicino è stata tagliata nel tetraedro vicino,
                //allora non si mette, si metterà in recupero conformità
                //Se non è stata tagliata, allora tra i vicini si mette il tetraedo (figlio o no) che contiene quella faccia
                //Per quanto riguarda 'second_face', si mettono i due vicini conformi
                GenericFace* face = meshPointer->Face(newCell_1.Face(3)->Id());
                //Con questo 'if' sai se il vicino è stato tagliato lungo questa faccia
                flag = 0;
                if (face->IsActive())
                {
                    unsigned int flag_2 = 0;
                    for (int position = 0; position < face->NumberOfCells() && flag_2 == 0; ++position)
                    {
                        if(face->Cell(position) != NULL)
                        {
                            unsigned int idTemp = face->Cell(position)->Id();
                            if (idTemp != newCell_1.Id() && idTemp != cell.Id())
                            {
                                GenericCell* cell_neig = meshPointer->Cell(idTemp);
                                while (cell_neig->HasChilds())
                                {
                                    unsigned int idChild = cell_neig->Child(0)->Id();
                                    GenericCell &child = *meshPointer->Cell(idChild);
                                    for (int faces = 0; faces < child.NumberOfFaces(); ++faces)
                                    {
                                        if (child.Face(faces)->Id() == face->Id())
                                        {
                                            ++flag;
                                            break;
                                        }
                                    }
                                    if (flag == 0)
                                    {
                                        idChild = cell_neig->Child(1)->Id();
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    else
                                    {
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    flag = 0;
                                }
                                /*flag = 0;
                                //Controllo per controllare se è già dentro a vicini (può succedere)
                                for (int faces = 0; faces < newCell_1.NumberOfFaces(); ++faces)
                                {
                                    unsigned int idTemp_2 = newCell_1.Face(faces)->Id();
                                    if (idTemp_2 == cell_neig.Id())
                                    {
                                        ++flag;
                                        break;
                                    }
                                }
                                if (flag == 0)
                                {*/
                                ++flag_2;
                                newCell_1.AddCell(cell_neig);
                                for (int cel = 0; cel < cell_neig->NumberOfCells(); cel++)
                                {
                                    if(cell_neig->Cell(cel) != NULL)
                                    {
                                        idTemp = cell_neig->Cell(cel)->Id();
                                        if (idTemp == cell.Id())
                                        {
                                            cell_neig->InsertCell(&newCell_1, cel);
                                        }
                                    }
                                }
                                //}
                                flag = 0;
                            }
                        }
                    }
                }
                face = meshPointer->Face(newCell_2.Face(3)->Id());
                if (face->IsActive())
                {
                    unsigned int flag_2 = 0;
                    for (int position = 0; position < face->NumberOfCells() && flag_2 == 0; ++position)
                    {
                        if(face->Cell(position) != NULL)
                        {
                            unsigned int idTemp = face->Cell(position)->Id();
                            if (idTemp != newCell_2.Id() && idTemp != cell.Id())
                            {
                                GenericCell* cell_neig = meshPointer->Cell(idTemp);
                                while (cell_neig->HasChilds())
                                {
                                    unsigned int idChild = cell_neig->Child(0)->Id();
                                    GenericCell &child = *meshPointer->Cell(idChild);
                                    for (int faces = 0; faces < child.NumberOfFaces(); ++faces)
                                    {
                                        if (child.Face(faces)->Id() == face->Id())
                                        {
                                            ++flag;
                                            break;
                                        }
                                    }
                                    if (flag == 0)
                                    {
                                        idChild = cell_neig->Child(1)->Id();
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    else
                                    {
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    flag = 0;
                                }
                                /*flag = 0;
                                //Controllo per controllare se è già dentro a vicini (può succedere)
                                for (int faces = 0; faces < newCell_2.NumberOfFaces(); ++faces)
                                {
                                    unsigned int idTemp_2 = newCell_2.Face(faces)->Id();
                                    if (idTemp_2 == cell_neig.Id())
                                    {
                                        ++flag;
                                        break;
                                    }
                                }
                                if (flag == 0)
                                {*/
                                ++flag_2;
                                newCell_2.AddCell(cell_neig);
                                for (int cel = 0; cel < cell_neig->NumberOfCells(); cel++)
                                {
                                    if(cell_neig->Cell(cel) != NULL)
                                    {
                                        idTemp = cell_neig->Cell(cel)->Id();
                                        if (idTemp == cell.Id())
                                        {
                                            cell_neig->InsertCell(&newCell_2, cel);
                                        }
                                    }
                                }
                                //}
                                flag = 0;
                            }
                        }
                    }
                }
                //Si è recuperata la conformità
                GenericCell* cell_neig = meshPointer->Cell(second_newEdge.Cell(0)->Id());
                newCell_1.AddCell(cell_neig);
                cell_neig->AddCell(&newCell_1);
                cell_neig = meshPointer->Cell(second_newEdge.Cell(1)->Id());
                newCell_2.AddCell(cell_neig);
                cell_neig->AddCell(&newCell_2);
                newCell_1.AddCell(&newCell_2);
                newCell_2.AddCell(&newCell_1);
            }
            else
                //****************************************************************************************************
                //Ultimo caso
            {
                //Si prendono i figli della prima faccia
                unsigned int id_first_nF = first_face.Child(0)->Id();
                GenericFace &first_newFace_1 = *(meshPointer->Face(id_first_nF));
                id_first_nF = first_face.Child(1)->Id();
                GenericFace &first_newFace_2 = *(meshPointer->Face(id_first_nF));

                //Si prende nuovo lato che taglia 'first_face'
                id_first_nF = first_newFace_1.Edge(2)->Id();
                GenericEdge &first_newEdge = *(meshPointer->Edge(id_first_nF));

                //****************************************************************************************************
                //Si prendono i figli della seconda faccia che condivide il lato più lungo
                //Si sa già qual è: 'second_face'

                unsigned id_second_nF = second_face.Child(0)->Id();
                GenericFace &second_newFace_1 = *(meshPointer->Face(id_second_nF));
                id_second_nF = second_face.Child(1)->Id();
                GenericFace &second_newFace_2 = *(meshPointer->Face(id_second_nF));

                //Si prende nuovo lato che taglia 'first_face'
                id_second_nF = second_newFace_1.Edge(2)->Id();
                GenericEdge &second_newEdge = *(meshPointer->Edge(id_second_nF));

                //****************************************************************************************************
                //Da qui si crea la terza nuova faccia

                GenericFace &third_newFace = *(meshPointer->CreateFace());
                meshPointer->AddFace(&third_newFace);
                third_newFace.AllocateCells(2);
                //Saranno inseriti i tetredri figli

                //Si aggiungono i punti alla faccia (si conoscono già) e viceversa
                third_newFace.InitializePoints(3);
                third_newFace.AddPoint(&newPoint);
                third_newFace.AddPoint(&first_opposite_point);
                third_newFace.AddPoint(&second_opposite_point);
                first_opposite_point.AddFace(&third_newFace);
                second_opposite_point.AddFace(&third_newFace);
                newPoint.AddFace(&third_newFace);

                //Si aggiungono i lati alla terza faccia e viceversa
                third_newFace.InitializeEdges(3);
                third_newFace.AddEdge(&first_newEdge);
                third_newFace.AddEdge(&second_newEdge);
                first_newEdge.AddFace(&third_newFace);
                second_newEdge.AddFace(&third_newFace);
                //Si deve andare a prendere il terzo lato che non è mai stato trovato (quello con 'first_opposite_point' e 'second_opposite_point'
                //come vertici)
                for (int position = 0; position < cell.NumberOfEdges(); ++position)
                {
                    unsigned int idTemp = cell.Edge(position)->Id();
                    GenericEdge &edge = *meshPointer->Edge(idTemp);
                    if (edge.Point(0)->Id() == id_first_opposite_point && edge.Point(1)->Id() == id_second_opposite_point)
                    {
                        third_newFace.AddEdge(&edge);
                        edge.AddFace(&third_newFace);
                        break;
                    }
                    else if (edge.Point(1)->Id() == id_first_opposite_point && edge.Point(0)->Id() == id_second_opposite_point)
                    {
                        third_newFace.AddEdge(&edge);
                        edge.AddFace(&third_newFace);
                        break;
                    }
                }

                //****************************************************************************************************
                //Da qui si creano i due tetraedri figli

                cell.SetState(false);
                cell.InitializeChilds(2);

                //Si crea il primo figlio e si settano alcuni suoi membri
                GenericCell &newCell_1 = *meshPointer->CreateCell();
                meshPointer->AddCell(&newCell_1);
                cell.AddChild(&newCell_1);
                newCell_1.SetFather(&cell);
                newCell_1.InheritPropertiesByFather();
                newCell_1.InitializeFaces(4);
                newCell_1.InitializeCells(4);
                newCell_1.InitializePoints(4);
                newCell_1.InitializeEdges(6);


                //Si crea il secondo figlio e si settano alcuni suoi membri
                GenericCell &newCell_2 = *meshPointer->CreateCell();
                meshPointer->AddCell(&newCell_2);
                cell.AddChild(&newCell_2);
                newCell_2.SetFather(&cell);
                newCell_2.InheritPropertiesByFather();
                newCell_2.InitializeFaces(4);
                newCell_2.InitializeCells(4);
                newCell_2.InitializePoints(4);
                newCell_2.InitializeEdges(6);


                //Si aggiungono la facce e viceversa
                newCell_1.AddFace(&third_newFace);
                newCell_2.AddFace(&third_newFace);
                third_newFace.AddCell(&newCell_1);
                third_newFace.AddCell(&newCell_2);
                //Si aggiungono le altre due facce coerentemente con il codice precedente (si sa già quali facce appartengono a quali figli:
                //si vuole che 'newCell_1' abbia 'newEdge_1', first_newFace_1', .....)
                newCell_1.AddFace(&first_newFace_1);
                newCell_2.AddFace(&first_newFace_2);
                first_newFace_1.AddCell(&newCell_1);
                first_newFace_2.AddCell(&newCell_2);
                newCell_1.AddFace(&second_newFace_1);
                newCell_2.AddFace(&second_newFace_2);
                second_newFace_1.AddCell(&newCell_1);
                second_newFace_2.AddCell(&newCell_2);
                flag = 0;
                //Si trova la quarta faccia per entrambi (facce mai toccate) e le si inserisce (e viceversa)
                for (unsigned int position = 0; position < cell.NumberOfFaces(); ++position)
                {
                    if (flag == 2)
                        break;
                    unsigned int idTemp = cell.Face(position)->Id();
                    if (idTemp != id_first_face && idTemp != id_second_face)
                    {
                        GenericFace &otherFace = *meshPointer->Face(idTemp);
                        for (int vertex = 0; vertex < otherFace.NumberOfPoints(); ++vertex)
                        {
                            unsigned int idTemp_2 = otherFace.Point(vertex)->Id();
                            if (idTemp_2 == id_first_point)
                            {
                                ++flag;
                                newCell_1.AddFace(&otherFace);
                                otherFace.AddCell(&newCell_1);
                                break;
                            }
                            else if (idTemp_2 == id_second_point)
                            {
                                ++flag;
                                newCell_2.AddFace(&otherFace);
                                otherFace.AddCell(&newCell_2);
                                break;
                            }
                        }
                    }
                }
                if(flag != 2)
                {
                    Output::PrintErrorMessage("The other two faces for the children of the cell with id %d were not found\n", false, cell.Id());
                    return Output::GenericError;
                }


                //Si aggiungono i lati ai tetraedri figli (e viceversa) coerentemente con il codice precedente
                newCell_1.AddEdge(&first_newEdge);
                newCell_1.AddEdge(&second_newEdge);
                first_newEdge.AddCell(&newCell_1);
                second_newEdge.AddCell(&newCell_1);
                newCell_2.AddEdge(&first_newEdge);
                newCell_2.AddEdge(&second_newEdge);
                first_newEdge.AddCell(&newCell_2);
                second_newEdge.AddCell(&newCell_2);
                newCell_1.AddEdge(&newEdge_1);
                newCell_2.AddEdge(&newEdge_2);
                newEdge_1.AddCell(&newCell_1);
                newEdge_2.AddCell(&newCell_2);
                //A partire da codice precedente
                GenericEdge* otherEdge = meshPointer->Edge(first_newFace_1.Edge(1)->Id());
                newCell_1.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_1);
                otherEdge = meshPointer->Edge(first_newFace_2.Edge(1)->Id());
                newCell_2.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_2);
                otherEdge = meshPointer->Edge(second_newFace_1.Edge(1)->Id());
                newCell_1.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_1);
                otherEdge = meshPointer->Edge(second_newFace_2.Edge(1)->Id());
                newCell_2.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_2);
                otherEdge = meshPointer->Edge(third_newFace.Edge(2)->Id());
                newCell_1.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_1);
                newCell_2.AddEdge(otherEdge);
                otherEdge->AddCell(&newCell_2);


                //Si aggiungono i punti ai tetraedri e viceversa
                newCell_1.AddPoint(&newPoint);
                newCell_2.AddPoint(&newPoint);
                newPoint.AddCell(&newCell_1);
                newPoint.AddCell(&newCell_2);
                newCell_1.AddPoint(&first_point);
                first_point.AddCell(&newCell_1);
                newCell_2.AddPoint(&second_point);
                second_point.AddCell(&newCell_2);
                newCell_1.AddPoint(&first_opposite_point);
                newCell_2.AddPoint(&first_opposite_point);
                first_opposite_point.AddCell(&newCell_1);
                first_opposite_point.AddCell(&newCell_2);
                newCell_1.AddPoint(&second_opposite_point);
                newCell_2.AddPoint(&second_opposite_point);
                second_opposite_point.AddCell(&newCell_1);
                second_opposite_point.AddCell(&newCell_2);

                //Si aggiungono i tetraedri vicini del padre (ancora conformi ai figli) ai figli coerentemente con la costruzione dei figli
                //Se la faccia (non con il lato più lungo della cella padre) che è in comune con il tetraedro vicino è stata tagliata nel tetraedro vicino,
                //allora non si mette, si metterà in recupero conformità
                //Se non è stata tagliata, allora tra i vicini si mette il tetraedo (figlio o no) che contiene quella faccia
                //Per quanto riguarda 'second_face' e 'first_face', si mettono i due vicini conformi
                GenericFace* face = meshPointer->Face(newCell_1.Face(3)->Id());
                //Con questo 'if' sai se il vicino è stato tagliato lungo questa faccia
                flag = 0;
                if (face->IsActive())
                {
                    unsigned int flag_2 = 0;
                    for (int position = 0; position < face->NumberOfCells() && flag_2 == 0; ++position)
                    {
                        if(face->Cell(position) != NULL)
                        {
                            unsigned int idTemp = face->Cell(position)->Id();
                            if (idTemp != newCell_1.Id() && idTemp != cell.Id())
                            {
                                GenericCell* cell_neig = meshPointer->Cell(idTemp);
                                while (cell_neig->HasChilds()) {
                                    unsigned int idChild = cell_neig->Child(0)->Id();
                                    GenericCell &child = *meshPointer->Cell(idChild);
                                    for (int faces = 0; faces < child.NumberOfFaces(); ++faces)
                                    {
                                        if (child.Face(faces)->Id() == face->Id())
                                        {
                                            ++flag;
                                            break;
                                        }
                                    }
                                    if (flag == 0)
                                    {
                                        idChild = cell_neig->Child(1)->Id();
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    else
                                    {
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    flag = 0;
                                }
                                /*flag = 0;
                                //Controllo per controllare se è già dentro a vicini (può succedere)
                                for (int faces = 0; faces < newCell_1.NumberOfFaces(); ++faces)
                                {
                                    unsigned int idTemp_2 = newCell_1.Face(faces)->Id();
                                    if (idTemp_2 == cell_neig.Id())
                                    {
                                        ++flag;
                                        break;
                                    }
                                }
                                if (flag == 0)
                                {*/
                                ++flag_2;
                                newCell_1.AddCell(cell_neig);
                                for (int cel = 0; cel < cell_neig->NumberOfCells(); cel++)
                                {
                                    if(cell_neig->Cell(cel) != NULL)
                                    {
                                        idTemp = cell_neig->Cell(cel)->Id();
                                        if (idTemp == cell.Id())
                                        {
                                            cell_neig->InsertCell(&newCell_1, cel);
                                        }
                                    }
                                }
                                //}
                                flag = 0;
                            }
                        }
                    }
                }
                face = meshPointer->Face(newCell_2.Face(3)->Id());
                if (face->IsActive())
                {
                    unsigned int flag_2 = 0;
                    for (int position = 0; position < face->NumberOfCells() && flag_2 == 0; ++position)
                    {
                        if(face->Cell(position) != NULL)
                        {
                            unsigned int idTemp = face->Cell(position)->Id();
                            if (idTemp != newCell_2.Id() && idTemp != cell.Id())
                            {
                                GenericCell* cell_neig = meshPointer->Cell(idTemp);
                                while (cell_neig->HasChilds())
                                {
                                    unsigned int idChild = cell_neig->Child(0)->Id();
                                    GenericCell &child = *meshPointer->Cell(idChild);
                                    for (int faces = 0; faces < child.NumberOfFaces(); ++faces)
                                    {
                                        if (child.Face(faces)->Id() == face->Id())
                                        {
                                            ++flag;
                                            break;
                                        }
                                    }
                                    if (flag == 0)
                                    {
                                        idChild = cell_neig->Child(1)->Id();
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    else
                                    {
                                        cell_neig = meshPointer->Cell(idChild);
                                    }
                                    flag = 0;
                                }
                                /*flag = 0;
                                //Controllo per controllare se è già dentro a vicini (può succedere)
                                for (int faces = 0; faces < newCell_2.NumberOfFaces(); ++faces)
                                {
                                    unsigned int idTemp_2 = newCell_2.Face(faces)->Id();
                                    if (idTemp_2 == cell_neig.Id()) {
                                        ++flag;
                                        break;
                                    }
                                }
                                if (flag == 0)
                                {*/
                                ++flag_2;
                                newCell_2.AddCell(cell_neig);
                                for (int cel = 0; cel < cell_neig->NumberOfCells(); cel++)
                                {
                                    if(cell_neig->Cell(cel) != NULL)
                                    {
                                        idTemp = cell_neig->Cell(cel)->Id();
                                        if (idTemp == cell.Id())
                                        {
                                            cell_neig->InsertCell(&newCell_2, cel);
                                        }
                                    }
                                }
                                //}
                                flag = 0;
                            }
                        }
                    }
                }
                //Si è recuperata la conformità
                GenericCell* cell_neig = meshPointer->Cell(second_newEdge.Cell(0)->Id());
                newCell_1.AddCell(cell_neig);
                cell_neig->AddCell(&newCell_1);
                cell_neig = meshPointer->Cell(second_newEdge.Cell(1)->Id());
                newCell_2.AddCell(cell_neig);
                cell_neig->AddCell(&newCell_2);
                cell_neig = meshPointer->Cell(first_newEdge.Cell(0)->Id());
                newCell_1.AddCell(cell_neig);
                cell_neig->AddCell(&newCell_1);
                cell_neig = meshPointer->Cell(first_newEdge.Cell(1)->Id());
                newCell_2.AddCell(cell_neig);
                cell_neig->AddCell(&newCell_2);
                newCell_1.AddCell(&newCell_2);
                newCell_2.AddCell(&newCell_1);
            }
        }

        return Output::Success;
    }


//****************************************************************************************************
//****************************************************************************************************
    const Output::ExitCodes RefinerTetra::AddIdCell(const unsigned int &idCell)
    {
        idCellToRefine.push_back(idCell);
        return Output::Success;
    }



//****************************************************************************************************
//****************************************************************************************************
    const Output::ExitCodes RefinerTetra::RecoverConformity()
    {
        //Se tetraedro è stato tagliato lungo lato, allora toglilo da vicini di quel lato se no viene ripreso: da guardare in cuttetra (così, ciò che devi scrivere sotto su figli: no caso di tetraedro inattivo
        //tagliato lungo lato che stai considernado qui, ha senso)

        //Si fa ciclo sui tetraedri che hanno lato (più lungo: lungo cui si taglia in 'CutTetra') in comune con un tetraedro precedentemente tagliato
        //passando attraverso i vicini del lato più lungo dei tetraedri tagliati precedentemente (ossia i tetraedri che condividono il lato)
        unsigned int idTemp;
        //unsigned int flag = 0;
        //int counter = 0;
        while (idEdgesCut.size() != 0)
        {
            //++counter;
            //printf("%d\n", counter);
            idTemp = idEdgesCut.front();
            GenericEdge &edge = *meshPointer->Edge(idTemp);
            for (int position = 0; position < edge.NumberOfCells(); ++position)
            {
                int idCell = edge.Cell(position)->Id();
                GenericCell &cell = *meshPointer->Cell(idCell);
                unsigned int idLongestEdge = 0;
                if (cell.IsActive())
                {
                    //printf("Cutting %d\n", cell.Id());
                    bool to_add;
                    CutTetra(cell, to_add);
                    FindLongestEdge(cell, idLongestEdge);
                    if(to_add == true)
                        idEdgesCut.push_back(idLongestEdge);
                    /*if (idLongestEdge != idTemp)
                    {
                        //list<unsigned int>::iterator it = find(idEdgesCut.begin(), idEdgesCut.end(), idLongestEdge);
                        //if(it == idEdgesCut.end())
                        if(!check_presence[idLongestEdge])
                        {
                            idEdgesCut.push_back(idLongestEdge);
                            check_presence[idLongestEdge] = true;
                        }*/
                        /*for (list<unsigned int>::iterator it = idEdgesCut.begin(); it != idEdgesCut.end(); ++it)
                        {
                            if (*it == idLongestEdge) {
                                flag = 1;
                                break;
                            }
                        }
                        if (flag == 0)
                        {
                            idEdgesCut.push_back(idLongestEdge);
                        }*/
                    //}
                }
                /*if(meshPointer->NumberOfCells() >= 348)
                {
                    printf("I punti sono:");
                    for(int i = 0; i < 4; ++i)
                    {
                        printf(" %d ", meshPointer->Cell(347)->Point(i)->Id());
                    }
                    printf("\n");
                    for(int i = 0; i < 4; ++i)
                    {
                        if(meshPointer->Cell(347)->Face(i)->Id() == 492)
                        {
                            printf("Yes ");
                            break;
                        }
                    }
                    printf("%d %d %d\n", meshPointer->Face(492)->Point(0)->Id(), meshPointer->Face(492)->Point(1)->Id(), meshPointer->Face(492)->Point(2)->Id());
                }*/
            }
            idEdgesCut.erase(idEdgesCut.begin());
        }

        CheckConformity();

        return Output::Success;
    }



//****************************************************************************************************
//****************************************************************************************************
    const Output::ExitCodes RefinerTetra::RefineMesh()
    {
        unsigned int idLongestEdge;
        while (idCellToRefine.size() != 0)
        {
            GenericCell &cellToRefine = *(meshPointer->Cell(idCellToRefine.front()));
            //printf("Cutting %d\n", cellToRefine.Id());
            bool to_add;
            RefinerTetra::CutTetra(cellToRefine, to_add);
            FindLongestEdge(cellToRefine, idLongestEdge);
            //unsigned int flag = 0;
            //list<unsigned int>::iterator it = find(idEdgesCut.begin(), idEdgesCut.end(), idLongestEdge);
            /*for (list<unsigned int>::iterator it = idEdgesCut.begin(); it != idEdgesCut.end(); ++it)
            {
                if (*it == idLongestEdge)
                {
                    flag = 1;
                    break;
                }
            }
            if (flag == 0)
            {
                idEdgesCut.push_back(idLongestEdge);
            }*/
            //if(it == idEdgesCut.end())
            if(to_add == true)
                idEdgesCut.push_back(idLongestEdge);
            /*if(!check_presence[idLongestEdge])
            {
                idEdgesCut.push_back(idLongestEdge);
                check_presence[idLongestEdge] = true;
            }*/
            idCellToRefine.erase(idCellToRefine.begin());
        }
        RefinerTetra::RecoverConformity();

        return Output::Success;
    }


//****************************************************************************************************
//****************************************************************************************************
    const Output::ExitCodes RefinerTetra::CheckConformity()
    {
        if(idEdgesCut.size() != 0)
        {
            Output::PrintErrorMessage("The are still tetrahedra to bring to conformity\n", false);
            return Output::GenericError;
        }

        unsigned int NumberOfFaces = meshPointer->NumberOfFaces();
        for(int position = 0; position < NumberOfFaces; ++position)
        {
            GenericFace& face = *meshPointer->Face(position);
            if(face.IsActive())
            {
                int count = 0;
                for(int cel = 0; cel < face.NumberOfCells(); ++cel)
                {
                    if(face.Cell(cel)!=NULL && face.Cell(cel)->IsActive())
                        ++count;
                }
                if(count > 2)
                {
                    Output::PrintErrorMessage("Conformity was not established: face with id %d belogns to more than two cells\n", false, face.Id());
                    return Output::GenericError;
                }
            }
        }

        unsigned int NumberOfCells = meshPointer->NumberOfCells();
        for(int cel = 0; cel < NumberOfCells; ++cel)
        {
            GenericCell& cell = *meshPointer->Cell(cel);
            if(cell.IsActive())
            {
                for(int edg = 0; edg < cell.NumberOfEdges(); ++edg)
                {
                    unsigned int idTemp = cell.Edge(edg)->Id();
                    if(!meshPointer->Edge(idTemp)->IsActive())
                    {
                        Output::PrintErrorMessage("Conformity was not established: edge with id %d in cell with id %d is not active\n", false, idTemp, cell.Id());
                        return Output::GenericError;
                    }
                }

                for(int fac = 0; fac < cell.NumberOfFaces(); ++fac)
                {
                    unsigned int idTemp = cell.Face(fac)->Id();
                    if(!meshPointer->Face(idTemp)->IsActive())
                    {
                        Output::PrintErrorMessage("Conformity was not established: face with id %d in cell with id %d is not active\n", false, idTemp, cell.Id());
                        return Output::GenericError;
                    }
                }
            }

        }

        return Output::Success;
    }
}
