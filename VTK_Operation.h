// VTK_Operation.h : interface of the VTK operation class
//
// Author: Jungho Park (jhpark16@gmail.com)
// Date: May 2015
// Description: 
//
/////////////////////////////////////////////////////////////////////////////

#pragma once

#include <vtkSmartPointer.h>
#include <vtkActor2D.h>
#include <vtkActor.h>
#include <vtkAxesActor.h>
#include <vtkTransform.h>
#include <vtkBMPReader.h>
#include <vtkCamera.h>
#include <vtkCaptionActor2D.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkConeSource.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkInteractorObserver.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkFloatArray.h>
#include <vtkHedgeHog.h>
#include <vtkHexahedron.h>
#include <vtkLine.h>
#include <vtkLookupTable.h>
#include <vtkMath.h>
#include <vtkNamedColors.h>
#include <vtkOpenGLRenderer.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPixel.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolygon.h>
#include <vtkPolyLine.h>
#include <vtkPolyVertex.h>
#include <vtkProperty.h>
#include <vtkPyramid.h>
#include "vtkScalarBarActor.h"
//#include "vtkScalarBarActor2.h"
//#include "vtkPVScalarBarActor.h"
//#include <vtkContext2DScalarBarActor.h>
#include <vtkStructuredGrid.h>
#include <vtkSphereSource.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkTextMapper.h>
#include <vtkTextProperty.h>
#include <vtkTexture.h>
#include <vtkTubeFilter.h>
#include <vtkQuad.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkTriangleStrip.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVoxel.h>
#include <vtkVertex.h>
#include <vtkWedge.h>
#include <vtkWin32OpenGLRenderWindow.h>
#include <vtkWin32RenderWindowInteractor.h>
#include <vtkXMLUnstructuredGridReader.h>

class VTK_Operation
{
public:

  VTK_Operation(HWND hwnd, HWND hwndParent);
  ~VTK_Operation();

  int CreateVTKObjects(HWND hwnd, HWND hwndParent);
  int DestroyVTKObjects();
  int StepObjects(int count);
  void WidgetInteractiveOff(void);
  void Render();

  void OnSize(CSize size);
  void OnTimer(HWND hWnd, UINT uMsg);
  void OnLButtonDblClk(HWND hWnd, UINT uMsg, CPoint point);
  void OnLButtonDown(HWND hWnd, UINT uMsg, CPoint point);
  void OnMButtonDown(HWND hWnd, UINT uMsg, CPoint point);
  void OnRButtonDown(HWND hWnd, UINT uMsg, CPoint point);
  void OnLButtonUp(HWND hWnd, UINT uMsg, CPoint point);
  void OnMButtonUp(HWND hWnd, UINT uMsg, CPoint point);
  void OnRButtonUp(HWND hWnd, UINT uMsg, CPoint point);
  void OnMouseMove(HWND hWnd, UINT uMsg, CPoint point);
  void OnChar(HWND hWnd, UINT nChar, UINT nRepCnt, UINT nFlags);
  void FileNew(void);
  void FileOpen(PHT3D_Model& mPHT3DM, CString fileName);
  void clear(void);

  //int DeleteObjects(void);

private:
  // renderer, window, and interactor
  vtkSmartPointer<vtkWin32OpenGLRenderWindow> renderWindow;
  vtkSmartPointer<vtkOpenGLRenderer> renderer;
  vtkSmartPointer<vtkWin32RenderWindowInteractor> interactor;
  // actor and mapper for text
  vtkSmartPointer<vtkActor2D> textActor2D;
  vtkSmartPointer<vtkTextMapper> textMapper;

  // actor and mapper for the model
  vtkSmartPointer<vtkActor> actor;
  //vtkSmartPointer<vtkDataSetMapper> mapper;

  // Axis
  vtkSmartPointer<vtkAxesActor> axes;

  // colormap
  vtkSmartPointer<vtkScalarBarActor> scalarBar;

  // Model reader, widget for axes triad, colors, and property
  
  vtkSmartPointer<vtkOrientationMarkerWidget> widget;
  CSize size;
};

