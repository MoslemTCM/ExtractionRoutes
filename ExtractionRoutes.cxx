/*========================================================================= */

// Équipe Télédetection pour les catastrophes majeures (TCM).

// Programme : Extraction des routes à partir des images satellitaires à THRS.

// Auteur : Moslem Ouled Sghaier

// Version : 0

/*========================================================================= */

#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "otbSFSTexturesImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkFlatStructuringElement.h"

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "PyramidTree.h"
#include <iostream>
#include "itkVector.h"
#include "itkBresenhamLine.h"
#include <vector>
#include <math.h>
#include <ctime>

#include "otbWrapperApplication.h" 
#include "otbWrapperApplicationRegistry.h"
#include "otbWrapperApplicationFactory.h"
#include "otbWrapperTags.h"

// Copier une image dans une autre
#include "itkPasteImageFilter.h"

//Utils
#include "itksys/SystemTools.hxx"
#include "itkListSample.h"

//otb filter
#include "otbExtractROI.h"

// Elevation handler
#include "otbWrapperElevationParametersHandler.h"

//using namespace std;

namespace otb
{
namespace Wrapper
{

class ExtractionRoutes : public Application
{

public:
	typedef ExtractionRoutes Self;
    typedef Application              Superclass;
    typedef itk::SmartPointer<Self>       Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;


#define PI 3.14159265

//using namespace std;
/////////////////////////////////////////////////////////////////////// déclaration des structures globales

 struct Point
{
int x;
int y;

 };


 struct Dic1
{
	int numero;
	int flag;
    std::vector<Point> p; 
};


 struct Dic
{
    int Level;
	std::vector<Dic1> D;

};


std::vector<Dic> Base; 

std::vector<Dic> Base1;

//////////////////////////////////////////////////////////////////////


typedef unsigned char CharPixelType; // IO

typedef itk::Image<CharPixelType, Dimension> ImageType;
//typedef otb::ImageFileReader<ImageType> ReaderType;
//typedef otb::ImageFileWriter<ImageType> WriterType;
typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;

////////////////////////////////////////////////////////////////////

int *Determiner(int x , int y)

{
  	int tab[2]; int x1=x/512; int y1=y/512;

	tab[0]=x1; tab[1]=y1;
    
 return tab; 
}

/** Standard macro */
    itkNewMacro(Self);
    itkTypeMacro(ExtractionRoutes, otb::Application);


// Les methodes ajoutées pour compiler le Qgis 
private:

void DoInit()
{

SetName("ExtractionRoutes"); // Nécessaire
SetDocName("ExtractionRoutes");
SetDocLongDescription("Un simple module pour l'extraction des routes");
SetDocLimitations("Les autres paramètres seront ajoutés plus tard");
SetDocAuthors("Moslem Ouled Sghaier");


AddParameter(ParameterType_InputImage,"in", "Input Image");
SetParameterDescription("in", "The input image");
AddParameter(ParameterType_OutputImage,"out", "Output Image");
SetParameterDescription("out","The output image");

}

void DoUpdateParameters()
{
	// Nothing to do here : all parameters are independent
}


void DoExecute()
{
  
	printf("Moslem va se  marrier CRIM");

  //const char * inputFilename = argv[1];
  //const char * outputFilename = argv[2];
  float variance = 45;
  float lowerThreshold = 1.5;
  float upperThreshold = 5.0;

  typedef unsigned char CharPixelType; // IO
  typedef double RealPixelType; // Operations
  typedef float FloatPixelType; // Operations
  typedef  UInt32ImageType MoslemImageType; // Operations
  const unsigned int Dimension = 2;

  typedef otb::Image<CharPixelType, Dimension> CharImageType;
  typedef otb::Image<RealPixelType, Dimension> RealImageType;
  typedef otb::Image<FloatPixelType, 2> FloatImageType;

  typedef otb::ImageFileReader<RealImageType>  ReaderType; 
  typedef otb::ImageFileWriter<RealImageType>  WriterType;

  typedef otb::SFSTexturesImageFilter<CharImageType, CharImageType> SFSFilterType;

  // Debut SFS-SD texture feature

  SFSFilterType::Pointer filter   = SFSFilterType::New(); 
 
  //ReaderType::Pointer    reader1   = ReaderType::New(); // je dois l'abondonner
  WriterType::Pointer    writerSD     = WriterType::New();

  typedef UInt32ImageType LabelImageType;

  typedef itk::RescaleIntensityImageFilter<UInt32ImageType,CharImageType> RescaleFilter0;
  RescaleFilter0::Pointer rescale0 = RescaleFilter0::New();
  rescale0->SetInput(GetParameterUInt32Image("in"));
  rescale0->UpdateLargestPossibleRegion();
  CharImageType::Pointer reader1=rescale0->GetOutput();
  //reader1->SetFileName(inputFilename);

  // Je voudrais récupirer les paramètres sur la longueur et la largeur de l'image

  reader1->UpdateOutputInformation();

  reader1->UpdateOutputData();

  CharImageType::SizeType size = reader1->GetLargestPossibleRegion().GetSize();

  filter->SetSpectralThreshold(25);
  filter->SetSpatialThreshold(100); 
  filter->SetNumberOfDirections(20); 
  filter->SetRatioMaxConsiderationNumber(5); 
  filter->SetAlpha(1.00);

  filter->SetInput(reader1);
  filter->Update();
  /// Fin de SFS SD

   // Debut Canny edge detection

  // premierement il faut caster le resultat en réel

  typedef otb::Image<unsigned char, 2> OutputImageType;
  typedef otb::ImageFileWriter<ImageType> OutputWriterType;

  OutputWriterType::Pointer outWriter = OutputWriterType::New();

  ////////////////////////////////////////////////////////////////////////ImageType::Pointer outputImage = ImageType::New();

  //outWriter->SetInput(rescaler->GetOutput());

  /// Début application de la dilatation 

  typedef itk::FlatStructuringElement<2>  StructuringElementType;
  typedef itk::GrayscaleDilateImageFilter<CharImageType,RealImageType,StructuringElementType>  GrayscaleDilateImageFilterType;
  GrayscaleDilateImageFilterType::Pointer dilateFilter = GrayscaleDilateImageFilterType::New();

  dilateFilter->SetInput(filter->GetOutput());
  StructuringElementType::RadiusType elementRadius;
  elementRadius.Fill(2);
  StructuringElementType structuringElement = StructuringElementType::Box(elementRadius);
  dilateFilter->SetKernel(structuringElement);

  // Canny

  typedef itk::CannyEdgeDetectionImageFilter<RealImageType,RealImageType> CannyFilter;
  CannyFilter::Pointer cannyFilter = CannyFilter::New();
  cannyFilter->SetInput(dilateFilter->GetOutput());
  cannyFilter->SetVariance(variance);

   // Seuillage
  cannyFilter->SetUpperThreshold( lowerThreshold );
  cannyFilter->SetLowerThreshold( upperThreshold );

  cannyFilter->UpdateLargestPossibleRegion();

  typedef itk::RescaleIntensityImageFilter<RealImageType,CharImageType> RescaleFilter1;
  RescaleFilter1::Pointer rescale1 = RescaleFilter1::New();

  //The output of an edge filter is 0 or 1
  rescale1->SetOutputMinimum(0);
  rescale1->SetOutputMaximum(255);

  rescale1->SetInput(cannyFilter->GetOutput());

  rescale1->Update();

  /// Fin de Canny edge detection

  // Inverser les pixels de l'image

  CharImageType::IndexType pixelIndex; 

  for (int i=0; i< size[0] ;i++)
  {
    for (int j=0; j< size[1] ;j++)
  {
	  pixelIndex[0]=i;pixelIndex[1]=j;

	  CharImageType::PixelType pixelValue = rescale1->GetOutput()->GetPixel(pixelIndex);
	  rescale1->GetOutput()->SetPixel(pixelIndex,std::abs(255-pixelValue));
	   
  }}

  // La partie extraction des routes, la transformée en Beamlet

  //Setting the IO


	ReaderType::Pointer reader = ReaderType::New();
	WriterType::Pointer writer = WriterType::New();

	//Création de l'image sortante

	CharImageType::Pointer outputfinal;
	outputfinal=rescale1->GetOutput();

    //Iterator
	//ImageType::Pointer inputImage;

    //inputImage = rescale1->GetOutput();

	std::cout << "La taille de l'image est" <<  size[0]  << "et" << size[1] << std::endl;
    int *tab = Determiner(size[0],size[1]);
    std::cout << "Le tab x est" <<  tab[0] << "et" << tab[1] << std::endl;
    printf("Le chargement des données est ok et la taille de l'image est %2d et %2d",size[0],size[1]);

	 if ( tab[0] >= 1 && tab[1] >= 1 ) // tawa kif zedtha
       {
     printf("On a %d et  %d  fenêtres de taille 512",tab[0],tab[1]); 
     int count=0;

	 for (int t=0;t< tab[0];t++)
	   for (int g=0;g< tab[1];g++)

	 {{
    
    printf("La valeur de count est %2d",count);
	count++;

    typedef otb::ExtractROI<ImageType::PixelType, ImageType::PixelType> FilterType; 
	FilterType::Pointer filter = FilterType::New();
    ImageType::ConstPointer inputImage;

	filter->SetStartX(t*512); 
    filter->SetStartY(g*512); 
    filter->SetSizeX(512); 
    filter->SetSizeY(512);

	filter->SetInput(rescale1->GetOutput());
    filter->Update();
	inputImage=filter->GetOutput();

	//initialisation de l'image de sortie
	ImageType::Pointer outputImage = ImageType::New();
	outputImage->SetRegions(inputImage->GetRequestedRegion());
	outputImage->CopyInformation(inputImage);
	outputImage->Allocate();
	outputImage->FillBuffer(255);
	//initialisation du premier niveau
	ImageType::RegionType::IndexType regionStart;
	ImageType::RegionType::SizeType regionSize;
	//w
	regionStart[0] = 0;
	//h
	regionStart[1] = 0;
	//w
	regionSize[0] = inputImage->GetRequestedRegion().GetSize()[0];
	//h
	regionSize[1] = inputImage->GetRequestedRegion().GetSize()[1];
	FortChildsTree *root;
	root = root->newFortChildsTree(regionStart, regionSize);

	divideAllImage(inputImage, outputImage, regionSize, regionStart, root);

	ImageType::SizeType inputSize = outputImage->GetLargestPossibleRegion().GetSize();

    double resultat= std::log10((double)inputSize[0]) / std::log10(2.0);

	cout << inputSize[0] << "taille de l'image" << resultat  <<  "Le nombre de niveau"<< std::endl ;

	
	for ( int p=0; p< (int)resultat ; p++ ) 
	{
	//ecrit un niveau
	std::list<FortChildsTree*> levelList;
	root->getLevelElements(root,p, &levelList);
	struct Dic d; // une structure pour remplir des informations
	d.Level=p;

	for (std::list<FortChildsTree*>::iterator levelIt = levelList.begin();levelIt != levelList.end(); ++levelIt) {
		
	regionScan(inputImage, outputImage, (*levelIt)->regionStart,(*levelIt)->regionSize,p,&d);
	
	}
     Base.push_back(d);
	}

	
	for(int y=0; y<Base.size() ; y++)
	cout << "L'echelle est : "  << Base.at(y).Level  << std::endl ;

	
	/// Première apparition seulement  !!!

    Base1=Base;

	premiere();
	seuillage();
    //angle();
    score();
    affiche(outputImage);

	typedef itk::RescaleIntensityImageFilter<ImageType,FloatImageType> RescaleFilter2;
    RescaleFilter2::Pointer rescale2 = RescaleFilter2::New();
	rescale2->SetInput(outputImage);

    
	rescale2->Update();

	// Fin du traitement

	typedef itk::RescaleIntensityImageFilter<FloatImageType,CharImageType> RescaleFilter5;
    RescaleFilter5::Pointer rescale5 = RescaleFilter5::New();
    rescale5->SetInput(rescale2->GetOutput());
    rescale5->UpdateLargestPossibleRegion();
	rescale5->SetOutputMaximum(0);
	rescale5->SetOutputMaximum(255);
	rescale5->Update();
	
	CharImageType::Pointer temp;
	temp=rescale5->GetOutput();

	typedef itk::PasteImageFilter <ImageType, ImageType > PasteImageFilterType;
	PasteImageFilterType::Pointer pasteFilter = PasteImageFilterType::New ();

	ImageType::IndexType destinationIndex;
    destinationIndex[0] = t*512;
    destinationIndex[1] = g*512;

	pasteFilter->SetSourceImage(temp);
    pasteFilter->SetDestinationImage(outputfinal);
    pasteFilter->SetSourceRegion(temp->GetLargestPossibleRegion());
    pasteFilter->SetDestinationIndex(destinationIndex);
	pasteFilter->Update();

	   }}
	SetParameterOutputImage<CharImageType>("out",outputfinal);
	 }
    else 
    printf("La taille de l'image est très petite pour appliquer des divisions"); 
}

void regionScan(ImageType::ConstPointer inputImage,
		ImageType::Pointer outputImage,
		ImageType::RegionType::IndexType regionStart,
		ImageType::RegionType::SizeType regionSize, int level,struct Dic *d) {
	std::vector<ImageType::IndexType> pointList;

	ImageType::RegionType region;
	region.SetSize(regionSize);
	region.SetIndex(regionStart);
	
	//iterator
	IteratorType regionIt(outputImage, region);
	IteratorType bresenhamIt(outputImage, region);
	ImageType::IndexType idx = regionIt.GetIndex();
	ImageType::IndexType firstPoint;
	
	//compteurs
	float total = 0;
	float error = 0;
	int beam = 0;
	

	for (regionIt.GoToBegin(); !regionIt.IsAtEnd(); ++regionIt) {

		if (isOnBoard(regionStart, regionSize, regionIt.GetIndex())) {
			if (regionIt.Get() != 0) {
				//cadriage
				//regionIt.Set(128);
			}
			idx = regionIt.GetIndex();
			if (inputImage->GetPixel(idx) == 0) {
				firstPoint = idx;
				for (bresenhamIt = regionIt; !bresenhamIt.IsAtEnd();
						++bresenhamIt) {
					idx = bresenhamIt.GetIndex();
					if (isOnBoard(regionStart, regionSize,
							bresenhamIt.GetIndex())) {
						if (inputImage->GetPixel(idx) == 0
								&& idx != firstPoint) {
							pointList = bresenham(firstPoint, idx);
							for (unsigned int i = 0; i < pointList.size();
									i++) {
								if (inputImage->GetPixel(pointList[i]) != 0) {
									error++;
								}
								total++;
							}
							if ((error / total) <= 0.40) {

								struct Dic1 d1;

								for (unsigned int i = 0; i < pointList.size();
										i++) {

                                    // il touchera plus
									//outputImage->SetPixel(pointList[i], 0);


									////////// remplir la structure ///////
									struct Point p1;
									
									p1.x =pointList[i].GetElement(0) ;
									p1.y =pointList[i].GetElement(1) ;
									d1.p.push_back(p1);
									
								}

							d1.numero=beam;
							beam++;
							d->D.push_back(d1);
							
							}
							
							error = 0;
							total = 0;
						}
					}
					pointList.clear();
				}
			}
		}
	}



}

void divideAllImage(ImageType::ConstPointer inputImage,
		ImageType::Pointer outputImage,
		ImageType::RegionType::SizeType regionSize,
		ImageType::RegionType::IndexType regionStart, FortChildsTree *root) {

	ImageType::RegionType::IndexType regionStartTemp;
	ImageType::RegionType::SizeType regionSizeTemp;

	if ((regionSize[0] / 2 != 1 && regionSize[1] / 2 != 1)) {
		//UpperLeft
		regionSizeTemp[0] = regionSize[0] / 2;
		regionSizeTemp[1] = regionSize[1] / 2;
		regionStartTemp = regionStart;
		divideAllImage(inputImage, outputImage, regionSizeTemp, regionStartTemp,
				root->insertTreeNode(root, regionStartTemp, regionSizeTemp, 1));
		//UpperRight
		regionStartTemp[0] = regionStart[0] + regionSizeTemp[0];
		regionStartTemp[1] = regionStart[1];
		divideAllImage(inputImage, outputImage, regionSizeTemp, regionStartTemp,
				root->insertTreeNode(root, regionStartTemp, regionSizeTemp, 2));
		//LowerLeft
		regionStartTemp[0] = regionStart[0];
		regionStartTemp[1] = regionStart[1] + regionSizeTemp[1];
		divideAllImage(inputImage, outputImage, regionSizeTemp, regionStartTemp,
				root->insertTreeNode(root, regionStartTemp, regionSizeTemp, 3));
		//LowerRight
		regionStartTemp[0] = regionStart[0] + regionSizeTemp[1];
		regionStartTemp[1] = regionStart[1] + regionSizeTemp[1];
		divideAllImage(inputImage, outputImage, regionSizeTemp, regionStartTemp,
				root->insertTreeNode(root, regionStartTemp, regionSizeTemp, 4));
	}
}

std::vector<ImageType::IndexType> bresenham(ImageType::IndexType start,
		ImageType::IndexType end) {

	int x1 = start[0];
	int y1 = start[1];
	int x2 = end[0];
	int y2 = end[1];

	ImageType::IndexType point;
	std::vector<ImageType::IndexType> pointList;

	int delta_x(x2 - x1);
	// if x1 == x2, then it does not matter what we set here
	signed char const ix((delta_x > 0) - (delta_x < 0));
	delta_x = std::abs(delta_x) << 1;

	int delta_y(y2 - y1);
	// if y1 == y2, then it does not matter what we set here
	signed char const iy((delta_y > 0) - (delta_y < 0));
	delta_y = std::abs(delta_y) << 1;

	point[0] = x1;
	point[1] = y1;
	pointList.push_back(point);
	if (delta_x >= delta_y) {
		// error may go below zero
		int error(delta_y - (delta_x >> 1));

		while (x1 != x2) {
			if ((error >= 0) && (error || (ix > 0))) {
				error -= delta_x;
				y1 += iy;
			}
			// else do nothing

			error += delta_y;
			x1 += ix;

			point[0] = x1;
			point[1] = y1;
			pointList.push_back(point);

		}
	} else {
		// error may go below zero
		int error(delta_x - (delta_y >> 1));

		while (y1 != y2) {
			if ((error >= 0) && (error || (iy > 0))) {
				error -= delta_y;
				x1 += ix;
			}
			// else do nothing

			error += delta_x;
			y1 += iy;

			point[0] = x1;
			point[1] = y1;
			pointList.push_back(point);
		}
	}
	return pointList;
}

bool isOnBoard(ImageType::RegionType::IndexType regionStart,
		ImageType::RegionType::SizeType regionSize,
		ImageType::IndexType pixel) {

	if (pixel[0] == regionStart[0]) {
		return true;
	}
	if (pixel[1] == regionStart[1]) {
		return true;
	}
	if ((pixel[0] + 1) % regionSize[0] == 0) {
		return true;
	}
	if ((pixel[1] + 1) % regionSize[1] == 0) {
		return true;
	}
	return false;
}

/////////// recherche d'apparition 

int exist(int i, int j, int x , int y)
{
    for (int q= 0 ;q<Base.size() && q< i ; q++ )

	{
		for (int w = 0 ;w< Base.at(q).D.size() && w<j ; w++ )

		{
		  
			Base.at(q).D.at(w);

			   for (int f = 0 ; f< Base.at(q).D.at(w).p.size()  ; f++ )

		 	   {
				  if (Base.at(q).D.at(w).p.at(f).x==x && Base.at(q).D.at(w).p.at(f).y==y )
                  return 1;
				

		 

		 	  }

       
		}



	}

return(0);

}

/////////// permet de determiner le nombre de beamlets à chaque étape


int compte()


{

int compteur=0;

for (int i= 0 ; i< Base.size(); i++ )

	{
		compteur+= Base.at(i).D.size();
    }   

return compteur;

}



////////// la première apparition 


void premiere ()

{

	int e=0;
	int e1=0;
	int supprime=0;


	for (int i= 0 ; i< Base.size(); i++ )

	{
		for (int j = 0 ; j< Base.at(i).D.size() ; j++ )

		{   
		    e1++;
			supprime=0;
			for (int k = 0 ; supprime!=1 && k< Base.at(i).D.at(j).p.size() ; k++ )

			{


				if (!Base.at(i).D.empty() && !Base.at(i).D.at(j).p.empty() )

				{

                 

				    if (  exist( i,j,Base.at(i).D.at(j).p.at(k).x,Base.at(i).D.at(j).p.at(k).y ) )

				    {
                    e++;
					Base.at(i).D.erase(Base.at(i).D.begin()+j);
					j--;
				    supprime=1;
				   }
                     
				}
				 
			}



		}
        

	}

	cout << "J'ai supprimé : " << e << std::endl;
	cout << "et ils sont  : " << e1 << std::endl;

}


void elimine_redondance()

{

  for (int i= 0 ; i< Base.size(); i++ )

	{
		for (int j = 0 ; j< Base.at(i).D.size() ; j++ )

		{
		   
			for (int k = 0 ; k< Base.at(i).D.at(j).p.size()-1 ; k++ )

			{
			
				///////////////////////////////////////////////////////////

				for (int k1 = k+1 ; k1< Base.at(i).D.at(j).p.size() ; k1++ )

			        {

						if (  (Base.at(i).D.at(j).p.at(k).x==Base.at(i).D.at(j).p.at(k1).x ) && (Base.at(i).D.at(j).p.at(k).y==Base.at(i).D.at(j).p.at(k1).y )  ){
							Base.at(i).D.at(j).p.erase(Base.at(i).D.at(j).p.begin()+k1);
					         k1--;
					    }
                    }

				///////////////////////////////////////////////////////////

			}

		}


  }




}




//////////// Ca permet de créer l'image résulante ///////////


void affiche(ImageType::Pointer outputImage)

{

	elimine_redondance();

	ImageType::IndexType point;
	int	e1=0;

	for (int i= 0 ; i< Base.size(); i++ )

	{
		for (int j = 0 ; j< Base.at(i).D.size() ; j++ )

		{
		    e1++;
			for (int k = 0 ; k< Base.at(i).D.at(j).p.size() ; k++ )

			{
             
			if (Base.at(i).D.at(j).p.size()>100){
            point[0]= Base.at(i).D.at(j).p.at(k).x;
            point[1]= Base.at(i).D.at(j).p.at(k).y;   
            outputImage->SetPixel( point , 0);
				}

			}

		}

	}

	cout << "et ils sont lors de l'affichage : " << e1 << std::endl;

}


void seuillage ()

{

	for (int i= 0 ; i< Base.size(); i++ )

	{
		for (int j = 0 ; j< Base.at(i).D.size() ; j++ )

		{
			 
			if (Base.at(i).D.at(j).p.size()<10)
			{
			 Base.at(i).D.erase(Base.at(i).D.begin()+j);
			  j--;
			}

		}


	}




}


Point verifier (int i,int j,Point p,Point p1)

{
int compteur =0;
//

 for (int k= i ; k< Base.size() ; k++ )

	{
		for (int h = j+1 ; h< Base.at(k).D.size() ; h++ )

		{


			int u =  Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).x - Base.at(i).D.at(j).p.at(0).x;
            int v =  Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).y - Base.at(i).D.at(j).p.at(0).y;
			int u1 = Base.at(k).D.at(h).p.at(Base.at(k).D.at(h).p.size()-1).x - Base.at(k).D.at(h).p.at(0).x;
            int v1 = Base.at(k).D.at(h).p.at(Base.at(k).D.at(h).p.size()-1).y - Base.at(k).D.at(h).p.at(0).y;

			int haut=(u*u1)+(v*v1) ;
			int bas=(sqrt((double)(u*u)+(v*v)))*(sqrt((double)(u1*u1)+(v1*v1)));


			int angle1 = acos((double)(haut)/(bas))* 180.0 / PI; ;

			int d= std::min (sqrt((double)((Base.at(i).D.at(j).p.at(0).x-Base.at(k).D.at(h).p.at(Base.at(k).D.at(h).p.size()-1).x)*(Base.at(i).D.at(j).p.at(0).x-Base.at(k).D.at(h).p.at(Base.at(k).D.at(h).p.size()-1).x))+
				((Base.at(i).D.at(j).p.at(0).y-Base.at(k).D.at(h).p.at(Base.at(k).D.at(h).p.size()-1).y)*(Base.at(i).D.at(j).p.at(0).y-Base.at(k).D.at(h).p.at(Base.at(k).D.at(h).p.size()-1).y))),
				sqrt((double)((Base.at(i).D.at(j).p.at(0).x-Base.at(k).D.at(h).p.at(0).x)*(Base.at(i).D.at(j).p.at(0).x-Base.at(k).D.at(h).p.at(0).x))+((Base.at(i).D.at(j).p.at(0).y-Base.at(k).D.at(h).p.at(0).y)*
				(Base.at(i).D.at(j).p.at(0).y-Base.at(k).D.at(h).p.at(0).y))));

			 int d1=std::min (sqrt((double)((Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).x-Base.at(k).D.at(h).p.at(Base.at(k).D.at(h).p.size()-1).x)*(Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).x-Base.at(k).D.at(h).p.at(Base.at(k).D.at(h).p.size()-1).x))+
				 ((Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).y-Base.at(k).D.at(h).p.at(Base.at(k).D.at(h).p.size()-1).y)*(Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).y-Base.at(k).D.at(h).p.at(Base.at(k).D.at(h).p.size()-1).y))),
				 sqrt((double)((Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).x-Base.at(k).D.at(h).p.at(0).x)*(Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).x-Base.at(k).D.at(h).p.at(0).x))+
				 ((Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).y-Base.at(k).D.at(h).p.at(0).y)*(Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).y-Base.at(k).D.at(h).p.at(0).y))));

			 int distance1=std::min(d,d1);

			// && distance1>0

			if (angle1<10 && distance1<15 ) // l'angle était 20 avant
			{
			 p.x=k;
             p.y=h;
             return p;
			
			}

			//cout << "L'angle entre les deux est " << angle1 << std::endl;
			compteur ++;
			


		}


    }
 //cout << "le nombre d'operation est  " << compteur << std::endl;

 p.x=0;
 p.y=0;

 return p;

}


void supprimer ()

{

for (int i= 0 ; i< Base.size(); i++ )

	{
		for (int j = 0 ; j< Base.at(i).D.size() ; j++ )

		{
			 
			if (Base.at(i).D.at(j).flag != 1)
			{
			 Base.at(i).D.erase(Base.at(i).D.begin()+j);
			  j--;
			}
		}
	}
}




void angle()

{


for (int i= 0 ; i< Base.size(); i++ )

	{
		for (int j = 0 ; j< Base.at(i).D.size() ; j++ )

		{
          
          if (Base.at(i).D.at(j).flag!=1)
		  {
          
          struct Point p=verifier (i,j, Base.at(i).D.at(j).p.at(0),Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1)  );

          if (p.x!=0 && p.y!=0)
		  {

			  Base.at(i).D.at(j).flag=1;
			  Base.at(p.x).D.at(p.y).flag=1;
			  //cout << "Ca marche !!!" << std::endl;
		  }

		  }

		}


    }

supprimer();

}

int verifier1(int i,int j,int k,int h)
{
// vérifier le segment vérifie les conditions
	        
	        // les pixels du début 
            int u =  Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).x - Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-2).x;
            int v =  Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).y - Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-2).y;
			int u1 = Base1.at(k).D.at(h).p.at(Base1.at(k).D.at(h).p.size()-1).x - Base1.at(k).D.at(h).p.at(0).x;
            int v1 = Base1.at(k).D.at(h).p.at(Base1.at(k).D.at(h).p.size()-1).y - Base1.at(k).D.at(h).p.at(0).y;

			// les pixels du fin
			int u2 =  Base.at(i).D.at(j).p.at(1).x - Base.at(i).D.at(j).p.at(0).x;
            int v2 =  Base.at(i).D.at(j).p.at(1).y - Base.at(i).D.at(j).p.at(0).y;
			int u3 = Base1.at(k).D.at(h).p.at(Base1.at(k).D.at(h).p.size()-1).x - Base1.at(k).D.at(h).p.at(0).x;
            int v3 = Base1.at(k).D.at(h).p.at(Base1.at(k).D.at(h).p.size()-1).y - Base1.at(k).D.at(h).p.at(0).y;

			// formule 1
			int haut=(u*u1)+(v*v1) ;
			int bas=(sqrt((double)(u*u)+(v*v)))*(sqrt((double)(u1*u1)+(v1*v1)));

			// formule 2
			int haut1=(u2*u3)+(v2*v3) ;
			int bas1=(sqrt((double)(u2*u2)+(v2*v2)))*(sqrt((double)(u3*u3)+(v3*v3)));

			int angle1 = acos((double)(haut)/(bas))* 180.0 / PI;
			int angle2 = acos((double)(haut1)/(bas1))* 180.0 / PI;


			int d= std::min (sqrt((double)((Base.at(i).D.at(j).p.at(0).x-Base1.at(k).D.at(h).p.at(Base1.at(k).D.at(h).p.size()-1).x)*(Base.at(i).D.at(j).p.at(0).x-Base1.at(k).D.at(h).p.at(Base1.at(k).D.at(h).p.size()-1).x))+
				((Base.at(i).D.at(j).p.at(0).y-Base1.at(k).D.at(h).p.at(Base1.at(k).D.at(h).p.size()-1).y)*(Base.at(i).D.at(j).p.at(0).y-Base1.at(k).D.at(h).p.at(Base1.at(k).D.at(h).p.size()-1).y))),
				sqrt((double)((Base.at(i).D.at(j).p.at(0).x-Base1.at(k).D.at(h).p.at(0).x)*(Base.at(i).D.at(j).p.at(0).x-Base1.at(k).D.at(h).p.at(0).x))+((Base.at(i).D.at(j).p.at(0).y-Base1.at(k).D.at(h).p.at(0).y)*
				(Base.at(i).D.at(j).p.at(0).y-Base1.at(k).D.at(h).p.at(0).y))));

			int d1=std::min (sqrt((double)((Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).x-Base1.at(k).D.at(h).p.at(Base1.at(k).D.at(h).p.size()-1).x)*(Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).x-Base1.at(k).D.at(h).p.at(Base1.at(k).D.at(h).p.size()-1).x))+
				 ((Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).y-Base1.at(k).D.at(h).p.at(Base1.at(k).D.at(h).p.size()-1).y)*(Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).y-Base1.at(k).D.at(h).p.at(Base1.at(k).D.at(h).p.size()-1).y))),
				 sqrt((double)((Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).x-Base1.at(k).D.at(h).p.at(0).x)*(Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).x-Base1.at(k).D.at(h).p.at(0).x))+
				 ((Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).y-Base1.at(k).D.at(h).p.at(0).y)*(Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).y-Base1.at(k).D.at(h).p.at(0).y))));

			int distance1=std::min(d,d1);

			
			//&&  distance1<2 
			if ( angle1 <20 &&  distance1<3   ) 
			{
			
             return 1;
			
			}

			else if ( angle2 <20 &&  distance1<3  ) 
			{
			
             return 2;
			
			}

return 0;
}



int recherche (int i, int j)

{

	for (int k= 0 ; k< Base1.size(); k++ )

	{
		for (int l = 0 ; l< Base1.at(k).D.size() ; l++ )

		{   
			// si la condition est vérifié mettre ajour le vecteur et supprimer l'élément 


			// 
			if (verifier1(i,j,k,l)==1)
			{
			
            // utiliser Breshenham pour relier 
             
            ImageType::IndexType deb;
			ImageType::IndexType fin;
             
			deb[0]=Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).x;
			deb[1]=Base.at(i).D.at(j).p.at(Base.at(i).D.at(j).p.size()-1).y;
			fin[0]=Base1.at(k).D.at(l).p.at(0).x;
			fin[1]=Base1.at(k).D.at(l).p.at(0).y;

			std::vector<ImageType::IndexType> vec; 
			vec=bresenham(deb,fin);

			if (vec.size() < 5)
			{
			for (unsigned int a=1 ; a < vec.size()-1 ; a++)
			{
			Point m;
			m.x=vec[a].GetElement(0);
			m.y=vec[a].GetElement(1);
			Base.at(i).D.at(j).p.push_back(m);
			}


            // Ajouter le nouveau 

			for (int d=0 ; d < Base1.at(k).D.at(l).p.size() ; d++  )
			{
            
			Base.at(i).D.at(j).p.push_back(Base1.at(k).D.at(l).p.at(d));
			}
            // Supprimer l'ancien

			Base1.at(k).D.erase(Base1.at(k).D.begin()+l);

			return 1;
			} }		 
            
			else if (verifier1(i,j,k,l)==2)
			{
			
            // utiliser Breshenham pour relier 
             
            ImageType::IndexType deb;
			ImageType::IndexType fin;
             
			deb[0]=Base.at(i).D.at(j).p.at(0).x;
			deb[1]=Base.at(i).D.at(j).p.at(0).y;
			fin[0]=Base1.at(k).D.at(l).p.at(Base1.at(k).D.at(l).p.size()-1).x;
			fin[1]=Base1.at(k).D.at(l).p.at(Base1.at(k).D.at(l).p.size()-1).y;

			std::vector<ImageType::IndexType> vec; 
			vec=bresenham(deb,fin);

			if (vec.size() < 3)
			{
			for (unsigned int a=1 ; a < vec.size()-1 ; a++)
			{
			Point m;
			m.x=vec[a].GetElement(0);
			m.y=vec[a].GetElement(1);
			Base.at(i).D.at(j).p.insert(Base.at(i).D.at(j).p.begin(),m);
			}


            // Ajouter le nouveau 

			for (int d=0 ; d < Base1.at(k).D.at(l).p.size() ; d++  )
			{
            
			Base.at(i).D.at(j).p.insert(Base.at(i).D.at(j).p.begin(),Base1.at(k).D.at(l).p.at(d));
			}
            // Supprimer l'ancien

			Base1.at(k).D.erase(Base1.at(k).D.begin()+l);

			return 1;
			}}		 
            


		}


	}


return 0;

}


void score()

{

	for (int i= 0 ; i< Base.size(); i++ )

	{
		for (int j = 0 ; j< Base.at(i).D.size() ; j++ )

		{
			 
			while(recherche(i,j))
			{}


		}

	}


}


       };
}
   
}
OTB_APPLICATION_EXPORT(otb::Wrapper::ExtractionRoutes);