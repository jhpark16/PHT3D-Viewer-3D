// VTK_Operation.cpp : implmentation of VTK operations to visualize a 3D model
//
// Author: Jungho Park (jhpark16@gmail.com)
// Date: May 2015
// Description: 
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "resource.h"

VTK_Operation::VTK_Operation(HWND hwnd, HWND hwndParent)
{
  CreateVTKObjects(hwnd, hwndParent);
}

void VTK_Operation::clear(void)
{
  // vtkSmartPointers are deleted by assigning nullptr 
  // (or any other object) to the pointer
#ifdef OpenVR
  actor.~vtkNew();
  renderWindow.~vtkNew();
  renderer.~vtkNew();
  interactor.~vtkNew();
  cam.~vtkNew();
  // destroy text actor and mapper
  // Additional variables
  //reader = nullptr;
  //mapper = nullptr;
  //colors = nullptr;
#else
  renderWindow.~vtkNew();
  renderer.~vtkNew();
  interactor.~vtkNew();
  // destroy text actor and mapper
  textMapper.~vtkNew();
  textActor2D.~vtkNew();
  scalarBar.~vtkNew();
  // Additional variables
  //reader = nullptr;
  //mapper = nullptr;
  //colors = nullptr;
  actor.~vtkNew();
  //backFace = nullptr;
#endif
}

VTK_Operation::~VTK_Operation()
{
  clear();
}

int VTK_Operation::CreateVTKObjects(HWND hwnd, HWND hwndParent)
{
  // We create the basic parts of a pipeline and connect them

#ifdef OpenVR
/*  vtkNew<vtkActor> actor;
  vtkNew<vtkOpenVRRenderer> renderer;
  vtkNew<vtkOpenVRRenderWindow> renderWindow;
  vtkNew<vtkOpenVRRenderWindowInteractor> interactor;
  vtkNew<vtkOpenVRCamera> cam;*/
/*  renderer = vtkSmartPointer<vtkOpenVRRenderer>::New();
  renderWindow = vtkSmartPointer<vtkOpenVRRenderWindow>::New();
  interactor = vtkSmartPointer<vtkOpenVRRenderWindowInteractor>::New();
  cam = vtkSmartPointer<vtkOpenVRCamera>::New();*/
#else
  /*
  actor = vtkSmartPointer<vtkActor>::New();
  renderer = vtkSmartPointer<vtkOpenGLRenderer>::New();
  renderWindow = vtkSmartPointer<vtkWin32OpenGLRenderWindow>::New();
  interactor = vtkSmartPointer<vtkWin32RenderWindowInteractor>::New();
  */
  renderWindow->Register(NULL);
#endif

#ifdef OpenVR
  renderWindow->AddRenderer(renderer.Get());
  renderer->AddActor(actor.Get());
  interactor->SetRenderWindow(renderWindow.Get());
  renderer->SetActiveCamera(cam.Get());

  vtkNew<vtkConeSource> cone;
  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputConnection(cone->GetOutputPort());
  actor->SetMapper(mapper.Get());

  //Reset camera to show the model at the centre
  renderer->ResetCamera();

#else
  renderWindow->AddRenderer(renderer);
  interactor->SetInstallMessageProc(0);
  // setup the parent window
  renderWindow->SetWindowId(hwnd);
  renderWindow->SetParentId(hwndParent);

  interactor->SetRenderWindow(renderWindow);
  interactor->Initialize();
  vtkNew<vtkInteractorStyleTrackballCamera> style;
  interactor->SetInteractorStyle(style);

  // Setup Text Actor
  //textMapper = vtkSmartPointer<vtkTextMapper>::New();
  textMapper->SetInput("MODFLOW Model");
  //textActor2D = vtkSmartPointer<vtkActor2D>::New();
  textActor2D->SetDisplayPosition(450, 550);
  textMapper->GetTextProperty()->SetFontSize(30);
  textMapper->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
  textActor2D->SetMapper(textMapper);
  // Add Axis
  //axes = vtkSmartPointer<vtkAxesActor>::New();
  axes->SetNormalizedLabelPosition(1.2, 1.2, 1.2);
  axes->GetXAxisCaptionActor2D()->SetHeight(0.025);
  axes->GetYAxisCaptionActor2D()->SetHeight(0.025);
  axes->GetZAxisCaptionActor2D()->SetHeight(0.025);
#ifndef OpenVR
  //widget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
  widget->SetOutlineColor(0.9300, 0.5700, 0.1300);
  widget->SetOrientationMarker(axes);
  widget->SetInteractor(interactor);
  widget->SetEnabled(1);
#endif
  // Prevent the interaction to fix the location of the coordinate axes triad
  // at the left bottom
  //widget->InteractiveOff(); // Trigger WM_SIZE
  widget->SetViewport(-0.8, -0.8, 0.25, 0.25);

  //Add the Actors
  renderer->AddActor(textActor2D);

  //renderer->SetBackground(colors->GetColor3d("Wheat").GetData());
  renderer->SetBackground(.1, .2, .4);

  //Reset camera to show the model at the centre
  renderer->ResetCamera();
#endif

  //renderer->SetBackground(0.0, 0.0, 0.25);
  return(1);
}


int VTK_Operation::DestroyVTKObjects()
{
  ATLTRACE("DestroyVTKObjects\n");
  return(1);
}

int VTK_Operation::StepObjects(int count)
{
  return(1);
}

void VTK_Operation::WidgetInteractiveOff(void)
{
#ifndef OpenVR
  widget->InteractiveOff(); // Trigger WM_SIZE
#endif
}

void VTK_Operation::Render()
{
  if (renderWindow) {
    renderWindow->Render();
    //interactor->Render();
  }
}

void VTK_Operation::OnSize(CSize size)
{
  if (renderWindow) {
    //int *val1;
    //val1 = renderWindow->GetSize();
    renderWindow->SetSize(size.cx, size.cy);
    if (textMapper && renderer->GetVTKWindow()) {
      int width = textMapper->GetWidth(renderer);
      this->size = size;
      textActor2D->SetDisplayPosition((size.cx - width) / 2, size.cy - 50);
      //interactor->UpdateSize(size.cx, size.cy);
      //interactor->SetSize(size.cx, size.cy);
    }
  }
}

void VTK_Operation::OnTimer(HWND hWnd, UINT uMsg)
{
#ifndef OpenVR
  interactor->OnTimer(hWnd, uMsg);
#endif
}

void VTK_Operation::OnLButtonDblClk(HWND hWnd, UINT uMsg, CPoint point)
{
#ifndef OpenVR
  interactor->OnLButtonDown(hWnd, uMsg, point.x, point.y, 1);
#endif
}

void VTK_Operation::OnLButtonDown(HWND hWnd, UINT uMsg, CPoint point)
{
#ifndef OpenVR
  interactor->OnLButtonDown(hWnd, uMsg, point.x, point.y, 0);
#endif
}

void VTK_Operation::OnMButtonDown(HWND hWnd, UINT uMsg, CPoint point)
{
#ifndef OpenVR
  interactor->OnMButtonDown(hWnd, uMsg, point.x, point.y, 0);
#endif
}

void VTK_Operation::OnRButtonDown(HWND hWnd, UINT uMsg, CPoint point)
{
#ifndef OpenVR
  interactor->OnRButtonDown(hWnd, uMsg, point.x, point.y, 0);
#endif
}

void VTK_Operation::OnLButtonUp(HWND hWnd, UINT uMsg, CPoint point)
{
#ifndef OpenVR
  interactor->OnLButtonUp(hWnd, uMsg, point.x, point.y);
#endif
}

void VTK_Operation::OnMButtonUp(HWND hWnd, UINT uMsg, CPoint point)
{
#ifndef OpenVR
  interactor->OnMButtonUp(hWnd, uMsg, point.x, point.y);
#endif
}

void VTK_Operation::OnRButtonUp(HWND hWnd, UINT uMsg, CPoint point)
{
#ifndef OpenVR
  interactor->OnRButtonUp(hWnd, uMsg, point.x, point.y);
#endif
}

void VTK_Operation::OnMouseMove(HWND hWnd, UINT uMsg, CPoint point)
{
#ifndef OpenVR
  interactor->OnMouseMove(hWnd, uMsg, point.x, point.y);
#endif
}

void VTK_Operation::OnChar(HWND hWnd, UINT nChar, UINT nRepCnt, UINT nFlags)
{
#ifndef OpenVR
  interactor->OnChar(hWnd, nChar, nRepCnt, nFlags);
#endif
}

void VTK_Operation::FileNew()
{
  renderer->RemoveActor(actor);
  renderer->RemoveActor(scalarBar);
}

void VTK_Operation::FileOpen(PHT3D_Model& mPHT3DM, CString fileName)
{
  // Model Data
  //read all the data from the file
  /*
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = 
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(_T("D:\\Study\\MODFLOW\\RegionalGroundwaterModelingwithMODFLOWandFlopy\\vtuFiles\\Model1_HD_Heads.vtu"));
  reader->Update();
  vtkSmartPointer<vtkUnstructuredGrid> output = reader->GetOutput();
  output->GetCellData()->SetScalars(output->GetCellData()->GetArray(0));
  */
  //vtkUnstructuredGrid construction
  vtkNew<vtkUnstructuredGrid> unstGrid;
  vtkNew<vtkFloatArray> fArr;
  vtkNew<vtkDoubleArray> dArr;
  vtkNew<vtkIntArray> iArr;
  vtkNew<vtkCellArray> cellArr;
  vtkNew<vtkIdTypeArray> cellIdType;
  vtkNew<vtkPoints> points;
  MODFLOWClass& mf = mPHT3DM.MODFLOW;
  double minVal=0, maxVal = 0;
  if (1) {
    OpenModel(mPHT3DM, fileName);
  int nCol = mf.DIS_NCOL + 1;
  int nRow = mf.DIS_NROW + 1;
  int nLay = mf.DIS_NLAY + 1;
  float *xLoc = new float[nCol]{};
  float *yLoc = new float[nRow]{};
  float *zLoc = new float[nRow*nCol]{};
  xLoc[0] = 0; yLoc[0] = 0;
  for (int i = 0; i < nCol-1; i++) {
    xLoc[i + 1] = xLoc[i] + mf.DIS_DELR[i];
  }
  for (int i = 0; i < nRow-1; i++) {
    yLoc[i + 1] = yLoc[i] + mf.DIS_DELC[i];
  }
  // Reverse Y for correct model display
  for (int i = 0; i < (nRow-1)/2; i++) {
    int iTmp = yLoc[i];
    yLoc[i] = yLoc[nRow - 1 - i];
    yLoc[nRow - 1 - i] = iTmp;
  }
  int nPoints = (mf.DIS_NLAY + 1)*(mf.DIS_NROW + 1)*(mf.DIS_NCOL + 1);
  unique_ptr<float[]> locations(new float[3 * nPoints]);
  //float *locations = new float[3 * nPoints]{};
  unique_ptr<int[]> connectivity(new int[9 * ((nLay - 1)*(nRow - 1)*(nCol - 1))]);
  unique_ptr<double[]> dVal(new double[((nLay - 1)*(nRow - 1)*(nCol - 1))]);
  
  //int *connectivity = new int[9*((nLay - 1)*(nRow - 1)*(nCol - 1))]{};
  //int *connectivity = new int[(nLay-1)*(nRow-1)*(nCol-1)] {};
  // The heights are not aligned at each nodal points. So, it is necessary to estimate the height
  // at each node using interpolation of heights
  CPPMatrix2<double> *elevation = &(mf.DIS_TOP);
  for (int k = 0; k < nLay; k++) {
    // Takes care of four corners
    if (k == 0) {
      // The elevation of the top layer is from the DIS_TOP
      elevation = &(mf.DIS_TOP);
    }
    else {
      // The bottom elevation is used for most layers
      elevation = &(mf.DIS_BOTMS[k - 1]);
    }
    zLoc[0 * nCol + 0] = (*elevation)[0][0];
    zLoc[0 * nCol + (nCol-1)] = (*elevation)[0][nCol-2];
    zLoc[(nRow - 1) * nCol + 0] = (*elevation)[nRow-2][0];
    zLoc[(nRow - 1) * nCol + (nCol-1)] = (*elevation)[nRow-2][nCol-2];
    // Interpolate edges along X direction 
    for (int i = 1; i < nCol-1; i++) {
      zLoc[0*nCol + i] = ((*elevation)[0][i-1] +
        (*elevation)[0][i])/2.0;
      zLoc[(nRow-1)*nCol + i] = ((*elevation)[nRow-2][i-1] +
        (*elevation)[nRow - 2][i])/2.0;
    }
    // Interpolate edges along Y direction 
    for (int j = 1; j < nRow - 1; j++) {
      zLoc[(j * nCol) + 0] = ((*elevation)[j-1][0] +
        (*elevation)[j][0]) / 2.0;
      zLoc[(j * nCol) + nCol-1] = ((*elevation)[j-1][nCol-2] +
        (*elevation)[j][nCol-2]) / 2.0;
    }
    // Interpolate the remaining part of the 3D model
    for (int j = 1; j < nRow-1; j++) {
      for (int i = 1; i < nCol-1; i++) {
        zLoc[(j * nCol) + i] = ((*elevation)[j - 1][i - 1] +
          (*elevation)[j][i - 1] + (*elevation)[j - 1][i] +
          (*elevation)[j][i]) / 4.0;
      }
    }
    // Set the X, Y, Z of all nodal points
    for (int j = 0; j < nRow; j++) {
      for (int i = 0; i < nCol; i++) {
        locations[3 * (k*(nRow*nCol) + j*(nCol)+i) + 0] = xLoc[i]; //X value
        locations[3 * (k*(nRow*nCol) + j*(nCol)+i) + 1] = yLoc[j]; //Y value
        locations[3 * (k*(nRow*nCol) + j*(nCol)+i) + 2] = zLoc[j*nCol+i]; //Z value
      }
    }
  }
  int numCells = 0;
  minVal = DBL_MAX;
  maxVal = -DBL_MAX;
  for (int k = 0; k < nLay-1; k++) {
    for (int j = 0; j < nRow-1; j++) {
      for (int i = 0; i < nCol-1; i++) {
        if (mf.BAS_IBOUND[k][j][i]) {
          connectivity[9 * numCells + 0] = 8;
          connectivity[9 * numCells + 1] = (k + 1)*(nRow*nCol) + (j + 1)*nCol + (i);
          connectivity[9 * numCells + 2] = (k + 1)*(nRow*nCol) + (j + 1)*nCol + (i + 1);
          connectivity[9 * numCells + 3] = (k+1)*(nRow*nCol) + (j)*nCol + (i+1);
          connectivity[9 * numCells + 4] = (k + 1)*(nRow*nCol) + (j)*nCol + (i);
          connectivity[9 * numCells + 5] = (k)*(nRow*nCol) + (j + 1)*nCol + (i);
          connectivity[9 * numCells + 6] = (k)*(nRow*nCol) + (j + 1)*nCol + (i + 1);
          connectivity[9 * numCells + 7] = (k)*(nRow*nCol) + (j)*nCol + (i + 1);
          connectivity[9 * numCells + 8] = (k)*(nRow*nCol) + (j)*nCol + (i);
          double tVal = mf.BAS_STRT[k][j][i];
          dVal[numCells] = tVal;
          if (tVal < minVal) minVal = tVal;
          if (tVal > maxVal) maxVal = tVal;
          numCells++;
        }
      }
    }
  }
  cellIdType->SetArray(connectivity.get(), numCells*9, 1);
  connectivity.release();
  fArr->SetNumberOfComponents(3); //3D data points
  fArr->SetArray(locations.get(), 3 * nPoints, 0);
  locations.release(); // Must release the pointer to prevent crash
  points->SetData(fArr); // The array should be allocated on the heap
  points->SetDataTypeToFloat();
    // Setup cell array using the raw data
    cellArr->SetCells(numCells, cellIdType);

    dArr->SetNumberOfComponents(1);
    dArr->SetArray(dVal.get(), numCells, 1);
    dVal.release(); // Must release the pointer to prevent crash
    for (int i = 1; i < numCells; i++) {
    }
    //dArr->SetTuple1(0, 0.1); dArr->SetTuple1(1, 0.5);
    // Assemble the unstructured grid
    unstGrid->SetPoints(points);
    unstGrid->SetCells(12, cellArr); // 12 is the cell type
    unstGrid->GetCellData()->SetScalars(dArr);
  }
  else {
    int nPoints = 12;
    int numCells = 2;
    //cellIdType->SetNumberOfValues(18);
    int *iCon = new int[numCells*9]{ 8, 0, 1, 2, 3, 4, 5, 6, 7,
      8, 4, 5, 6, 7, 8, 9, 10, 11 };
    cellIdType->SetArray(iCon, numCells * 9, 1);
    /*
    cellIdType->SetValue(0, 8);  cellIdType->SetValue(1, 0);  cellIdType->SetValue(2, 1);
    cellIdType->SetValue(3, 2);  cellIdType->SetValue(4, 3);  cellIdType->SetValue(5, 4);
    cellIdType->SetValue(6, 5);  cellIdType->SetValue(7, 6);  cellIdType->SetValue(8, 7);
    cellIdType->SetValue(9, 8);  cellIdType->SetValue(10, 4);  cellIdType->SetValue(11, 5);
    cellIdType->SetValue(12, 6);  cellIdType->SetValue(13, 7);  cellIdType->SetValue(14, 8);
    cellIdType->SetValue(15, 9);  cellIdType->SetValue(16, 10);  cellIdType->SetValue(17, 11);
    */
    cellArr->SetCells(numCells, cellIdType);
    // Set points
    double *(val[12]);
    float *fPos = new float[nPoints * 3]{ 0, 0, 0,  0, 0, 1,  0, 1, 1,  0, 1, 0,
      1, 0, 0,  1, 0, 1,  1, 1, 1,  1, 1, 0,
      2, 0, 0,  2, 0, 1,  2, 1, 1,  2, 1, 0 };
    fArr->SetNumberOfComponents(3); //3D data points
                                    //fArr->SetNumberOfTuples(12); // 12 set of 3D data points (not necessary)
    fArr->SetArray(fPos, nPoints*3, 0);
    points->SetData(fArr); // The array should be allocated on the heap
    points->SetDataTypeToFloat();
    /*
    for (int i = 0; i<12; i++)
    val[i] = points->GetPoint(i);
    points->SetNumberOfPoints(12);
    points->SetPoint(0, 0, 0, 0.);  points->SetPoint(1, 0, 0, 1.);  points->SetPoint(2, 0, 1, 1.);
    points->SetPoint(3, 0, 1, 0.);  points->SetPoint(4, 1, 0, 0.);  points->SetPoint(5, 1, 0, 1.);
    points->SetPoint(6, 1, 1, 1.);  points->SetPoint(7, 1, 1, 0.);  points->SetPoint(8, 2, 0, 0.);
    points->SetPoint(9, 2, 0, 1.);  points->SetPoint(10, 2, 1, 1.);  points->SetPoint(11, 2, 1, 0.);
    points->SetDataTypeToFloat();  // Point set should be float type
    for (int i = 0; i<12; i++)
    val[i] = points->GetPoint(i);
    */
    cellArr->SetCells(numCells, cellIdType);
    dArr->SetNumberOfComponents(1);
    //dArr->SetNumberOfTuples(2);
    minVal = 0, maxVal = 1;
    dArr->SetArray(new double[numCells] {0.1,0.5}, numCells, 1);
    //dArr->SetTuple1(0, 0.1); dArr->SetTuple1(1, 0.5);
    unstGrid->SetPoints(points);
    unstGrid->SetCells(12, cellArr); // 12 is the cell type
    unstGrid->GetCellData()->SetScalars(dArr);
  }
  // Setup Model Actor
  vtkNew<vtkDataSetMapper> mapper;
  //mapper->SetInputConnection(reader->GetOutputPort());
  mapper->SetInputData(unstGrid);
  //mapper->ScalarVisibilityOff();
  mapper->ScalarVisibilityOn();

  // Scale bar actor for the Colormap
  //scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
  scalarBar->SetLookupTable(mapper->GetLookupTable());
  scalarBar->SetUnconstrainedFontSize(true);
  scalarBar->SetTitle("Head (ft)");
  scalarBar->SetVerticalTitleSeparation(5);
  scalarBar->GetTitleTextProperty()->SetFontSize(20);
  scalarBar->GetTitleTextProperty()->ItalicOff();
  scalarBar->SetLabelFormat("%.1f");
  scalarBar->GetLabelTextProperty()->ItalicOff();
  scalarBar->GetLabelTextProperty()->SetFontSize(15);
  scalarBar->SetPosition(0.87, 0.2);
  scalarBar->SetHeight(0.5);
  scalarBar->SetWidth(0.11);
  scalarBar->SetNumberOfLabels(4);
  // Jet color scheme
  vtkNew<vtkLookupTable> jet;
  jet->SetNumberOfColors(257);
  jet->Build();
  for (int i = 0; i <= 64; i++) {
    jet->SetTableValue(i, 0, i / 64.0, 1, 1);
  }
  for (int i = 65; i <= 128; i++) {
    jet->SetTableValue(i, 0, 1, 1.0 - (i - 64.0) / 64.0, 1);
  }
  for (int i = 129; i <= 192; i++) {
    jet->SetTableValue(i, (i - 128.0) / 64.0, 1, 0, 1);
  }
  for (int i = 193; i <= 256; i++) {
    jet->SetTableValue(i, 1, 1 - (i - 192.0) / 64.0, 0, 1);
  }
  //jet->SetTableRange(3407.6445, 4673.021);
  mapper->SetLookupTable(jet);
  scalarBar->SetLookupTable(jet);

  mapper->SetColorModeToMapScalars();
  //mapper->SetScalarRange(3407.6445, 4673.021); // min max range cannot be switched
  mapper->SetScalarRange(minVal, maxVal); // min max range cannot be switched
  mapper->SetUseLookupTableScalarRange(false);
  // Set scalar mode
  mapper->SetScalarModeToUseCellData();
  //mapper->SetScalarModeToUsePointData();

  actor->SetMapper(mapper);
  // True for grid lines
  if (true) {
    actor->GetProperty()->SetLineWidth(0.001); // This should not be zero => OpenGL error (invalid value)
    actor->GetProperty()->EdgeVisibilityOn();
  }
  else
    actor->GetProperty()->EdgeVisibilityOff();

  // Backface setup
  vtkNew<vtkNamedColors> colors;
  vtkNew<vtkProperty> backFace;
  backFace->SetColor(colors->GetColor3d("tomato").GetData());
  actor->SetBackfaceProperty(backFace);

  renderer->AddActor(actor);
  renderer->AddActor(scalarBar);

  //Reset camera to show the model at the centre
  renderer->ResetCamera();
#ifdef OpenVR
  renderWindow->Render();
  interactor->Start();
#endif
}
