#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/freeglut.h>

#include <optixu/optixpp_namespace.h>
#include <optixu/optixu_math_stream_namespace.h>

#include <sutil.h>
#include <Arcball.h>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <math.h>
#include <stdint.h>

//#define TONE_MAPPING_POST_PROCESSING_STAGE

const char* const SAMPLE_NAME = "optixDLDenoiserVolumeRendering";
const int NUM_OF_NON_DENOISED_FRAMES = 5;

#include "voldatamanager.h"
VolumeDataManager voldatamanager;

int width = 768;
int height = 768;
bool use_pbo = true;
optix::Context context;

const char*  tutorial_ptx;

// Camera state
optix::float3 camera_up;
optix::float3 camera_lookat;
optix::float3 camera_eye;
optix::Matrix4x4 camera_rotate;
sutil::Arcball arcball;

// Mouse state
optix::int2 mouse_prev_pos;
int  mouse_button;

int show_buffer = 0;

optix::int3 vol_sizes;
optix::float3 vol_scaledsizes;

// denoiser variables
bool postprocessing_must_init = true;
bool camera_changed = true;
int frame_number = 1;

float denoise_blend = 0.f;

optix::CommandList commandListWithDenoiser;
optix::CommandList commandListWithoutDenoiser;
bool post_processing_stage_created = false;
#ifdef TONE_MAPPING_POST_PROCESSING_STAGE
optix::PostprocessingStage pps_tonemap;
#endif
optix::PostprocessingStage pps_denoiser;
optix::Buffer empty_buffer;
// will be empty, but we can use custom training data
optix::Buffer detraining_buffer;
optix::Buffer denoised_buffer;

optix::Buffer GetOutputBuffer ()
{
  return context["out_render_buffer"]->getBuffer();
}

optix::Buffer GetTonemappedBuffer ()
{
  return context["out_tonemap_buffer"]->getBuffer();
}

optix::Buffer GetAlbedoBuffer ()
{
  return context["out_albedo_buffer"]->getBuffer();
}

optix::Buffer GetNormalBuffer ()
{
  return context["out_normal_buffer"]->getBuffer();
}

optix::Buffer GetDenoisedBuffer ()
{
  return denoised_buffer;
}

void CreateOptixContext ()
{
  // OptiX Engine instance "Context"
  context = optix::Context::create();

  //context->setRayTypeCount(2);
  // Actually, we don't need 2 types of ray, we don't have shadow for tutorial0
  rtContextSetRayTypeCount(context->get(), 1);
  rtContextSetEntryPointCount(context->get(), 1);
  rtContextSetMaxTraceDepth(context->get(), 1);

  context["scene_epsilon"]->setFloat(1.e-4f);

  // render buffer
  optix::Buffer buffer_render = sutil::createInputOutputBuffer(context, RT_FORMAT_FLOAT4, width, height, use_pbo);
  context["out_render_buffer"]->set(buffer_render);
  // albedo buffer
  optix::Buffer buffer_albedo = sutil::createInputOutputBuffer(context, RT_FORMAT_FLOAT4, width, height, use_pbo);
  context["out_albedo_buffer"]->set(buffer_albedo);
  // normal buffer
  optix::Buffer buffer_normal = sutil::createInputOutputBuffer(context, RT_FORMAT_FLOAT4, width, height, use_pbo);
  context["out_normal_buffer"]->set(buffer_normal);

  // tonemapped buffer
  optix::Buffer buffer_tonemap = sutil::createInputOutputBuffer(context, RT_FORMAT_FLOAT4, width, height, use_pbo);
  context["out_tonemap_buffer"]->set(buffer_tonemap);

  // final denoised buffer
  denoised_buffer = sutil::createOutputBuffer(context, RT_FORMAT_FLOAT4, width, height, use_pbo);
  
  ///////////////////////////////////////////////
  // general empty buffer
  empty_buffer = context->createBuffer(RT_BUFFER_OUTPUT, RT_FORMAT_FLOAT4, 0, 0);
  // training data buffer
  detraining_buffer = context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_UNSIGNED_BYTE, 0);
  // use default value (without using custom training data)
  detraining_buffer->setSize(0);
  ///////////////////////////////////////////////

  // Ray generation program
  rtContextSetRayGenerationProgram(context->get(), 0, context->createProgramFromPTXString(tutorial_ptx, "RayGenerationProgram")->get());
  // Exception program
  rtContextSetExceptionProgram(context->get(), 0, context->createProgramFromPTXString(tutorial_ptx, "exception")->get());
  context["bad_color"]->setFloat(0.0f, 0.0f, 0.0f);
  // Miss program
  rtContextSetMissProgram(context->get(), 0, context->createProgramFromPTXString(tutorial_ptx, "miss")->get());
  context["bg_color"]->setFloat(optix::make_float3(0.0f, 0.0f, 0.0f));
}

void DestroyOptixContext ()
{
  if (context)
  {
    context->destroy();
    context = 0;
  }
}

void CreateSceneGeometry ()
{
  // Box geometry
  optix::Geometry box = context->createGeometry();
  box->setPrimitiveCount( 1u );
  const char *ptx = sutil::getPtxString( SAMPLE_NAME, "box.cu" );
  box->setBoundingBoxProgram(context->createProgramFromPTXString(ptx, "box_bounds"));
  box->setIntersectionProgram(context->createProgramFromPTXString(ptx, "box_intersect"));
  box["boxmin"]->setFloat(-vol_scaledsizes * 0.5f);
  box["boxmax"]->setFloat(vol_scaledsizes * 0.5f);

  // Box Material? just create but don't set anything relevant
  optix::Material box_matl = context->createMaterial();
  optix::Program box_ch = context->createProgramFromPTXString( tutorial_ptx, "closest_hit_radiance0");
  box_matl->setClosestHitProgram( 0, box_ch );

  // Create GIs for each piece of geometry
  std::vector<optix::GeometryInstance> gis;
  gis.push_back(context->createGeometryInstance( box, &box_matl, &box_matl+1 ) );

  // Place all in group
  optix::GeometryGroup geometrygroup = context->createGeometryGroup();
  geometrygroup->setChildCount( static_cast<unsigned int>(gis.size()) );
  geometrygroup->setChild( 0, gis[0] );
  
  geometrygroup->setAcceleration(context->createAcceleration("Trbvh") );
  //geometrygroup->setAcceleration( context->createAcceleration("NoAccel") );

  context["top_object"]->set( geometrygroup );
  context["top_shadower"]->set( geometrygroup );
}

void CreateTextureSamplers ()
{
  voldatamanager.ReadVolume("../../optixDLDenoiserVolumeRendering/data/Bucky.pvm");
  //voldatamanager.ReadVolume("../../optixDLDenoiserVolumeRendering/data/Engine.pvm");
  voldatamanager.ReadTransferFunction("../../optixDLDenoiserVolumeRendering/data/isovolume.tf1d");
  voldatamanager.SetVolumeTextureSamplerToOptix(context->get());
  voldatamanager.SetTransferFunctionTextureSamplerToOptix(context->get());
  voldatamanager.SetGradientTextureSamplerToOptix(context->get());

  vol_sizes = voldatamanager.GetVolumeSizes();
  vol_scaledsizes = voldatamanager.GetScaledBoundingBox();
}

void SetupPostProcessingStage ()
{
  // if post processing stage not created
  if (!post_processing_stage_created)
  {
#ifdef TONE_MAPPING_POST_PROCESSING_STAGE
    pps_tonemap = context->createBuiltinPostProcessingStage("TonemapperSimple");
#endif
    pps_denoiser = context->createBuiltinPostProcessingStage("DLDenoiser");
    if (detraining_buffer)
    {
      optix::Variable trainingbuff = pps_denoiser->declareVariable("training_data_buffer");
      trainingbuff->set(detraining_buffer);
    }
#ifdef TONE_MAPPING_POST_PROCESSING_STAGE
    pps_tonemap->declareVariable("input_buffer")->set(GetOutputBuffer());
    pps_tonemap->declareVariable("output_buffer")->set(GetTonemappedBuffer());
    pps_tonemap->declareVariable("exposure")->setFloat(2.0f);
    pps_tonemap->declareVariable("gamma")->setFloat(2.0f);

    pps_denoiser->declareVariable("input_buffer")->set(GetTonemappedBuffer());
#else 
    pps_denoiser->declareVariable("input_buffer")->set(GetOutputBuffer());
#endif

    pps_denoiser->declareVariable("output_buffer")->set(GetDenoisedBuffer());
    pps_denoiser->declareVariable("blend")->setFloat(denoise_blend);
    pps_denoiser->declareVariable("input_albedo_buffer");
    pps_denoiser->declareVariable("input_normal_buffer");

    optix::Variable albedoBuffer = pps_denoiser->queryVariable("input_albedo_buffer");
    albedoBuffer->set(GetAlbedoBuffer());
    optix::Variable normalBuffer = pps_denoiser->queryVariable("input_normal_buffer");
    normalBuffer->set(GetNormalBuffer());

    post_processing_stage_created = true;
  }

  if (commandListWithDenoiser)
  {
    commandListWithDenoiser->destroy();
    commandListWithoutDenoiser->destroy();
  }

  // Two command list will be created. One that contains 
  //  the denoiser and another that does not use the denoiser.
  commandListWithDenoiser = context->createCommandList();
  commandListWithDenoiser->appendLaunch(0, width, height);
#ifdef TONE_MAPPING_POST_PROCESSING_STAGE
  commandListWithDenoiser->appendPostprocessingStage(pps_tonemap, width, height);
#endif
  commandListWithDenoiser->appendPostprocessingStage(pps_denoiser, width, height);
  commandListWithDenoiser->finalize();

  commandListWithoutDenoiser = context->createCommandList();
  commandListWithoutDenoiser->appendLaunch(0, width, height);
#ifdef TONE_MAPPING_POST_PROCESSING_STAGE
  commandListWithoutDenoiser->appendPostprocessingStage(pps_tonemap, app.width, app.height);
#endif
  commandListWithoutDenoiser->finalize();

  postprocessing_must_init = false;
}

void UpdateCamera ()
{
  const float vfov = 60.0f;
  const float aspect_ratio = static_cast<float>(width) /
    static_cast<float>(height);

  optix::float3 camera_u, camera_v, camera_w;
  sutil::calculateCameraVariables(
    camera_eye, camera_lookat, camera_up, vfov, aspect_ratio,
    camera_u, camera_v, camera_w, true);

  const optix::Matrix4x4 frame = optix::Matrix4x4::fromBasis(
    optix::normalize(camera_u),
    optix::normalize(camera_v),
    optix::normalize(-camera_w),
    camera_lookat);
  const optix::Matrix4x4 frame_inv = frame.inverse();
  // Apply camera rotation twice to match old SDK behavior
  const optix::Matrix4x4 trans = frame * camera_rotate*camera_rotate*frame_inv;

  camera_eye = optix::make_float3(trans*optix::make_float4(camera_eye, 1.0f));
  camera_lookat = optix::make_float3(trans*optix::make_float4(camera_lookat, 1.0f));
  camera_up = optix::make_float3(trans*optix::make_float4(camera_up, 0.0f));

  sutil::calculateCameraVariables(
    camera_eye, camera_lookat, camera_up, vfov, aspect_ratio,
    camera_u, camera_v, camera_w, true);

  camera_rotate = optix::Matrix4x4::identity();

  if (camera_changed) // reset accumulation
    frame_number = 1;
  camera_changed = false;

  context["frame_number"]->setUint(frame_number);
  context["eye"]->setFloat(camera_eye);
  context["U"]->setFloat(camera_u);
  context["V"]->setFloat(camera_v);
  context["W"]->setFloat(camera_w);

  /////////////////////////////////////////////////////////////////////////
  // Store the inverse camera matrix to be multiplied by normals
  const optix::Matrix4x4 current_frame_inv = optix::Matrix4x4::fromBasis(
    normalize(camera_u),
    normalize(camera_v),
    normalize(-camera_w),
    camera_lookat).inverse();
  optix::Matrix3x3 normal_matrix = make_matrix3x3(current_frame_inv);

  context["normal_matrix"]->setMatrix3x3fv(false, normal_matrix.getData());
}

void glutDisplay()
{
  // Update camera and reset frame number if camera changed
  UpdateCamera();

  context["VolSizes"]->setFloat(vol_scaledsizes);

  if (postprocessing_must_init)
    SetupPostProcessingStage();

  // update blend between denoised image and the original image
  optix::Variable(pps_denoiser->queryVariable("blend"))->setFloat(denoise_blend);

  bool is_early_frame = (frame_number <= NUM_OF_NON_DENOISED_FRAMES);
  if (is_early_frame)
    commandListWithoutDenoiser->execute();
  else
    commandListWithDenoiser->execute();

  if (show_buffer == 0)
  {
    if (is_early_frame)
    {
#ifdef TONE_MAPPING_POST_PROCESSING_STAGE
      sutil::displayBufferGL(GetTonemappedBuffer(), BUFFER_PIXEL_FORMAT_DEFAULT, true);
#else
      sutil::displayBufferGL(GetOutputBuffer());
#endif
    }
    else
    {
      sutil::displayBufferGL(GetDenoisedBuffer());
    }
  }
  else if (show_buffer == 1)
    sutil::displayBufferGL(GetOutputBuffer());
  else if (show_buffer == 2)
    sutil::displayBufferGL(GetAlbedoBuffer());
  else if (show_buffer == 3)
    sutil::displayBufferGL(GetNormalBuffer());
  else if (show_buffer == 4)
    sutil::displayBufferGL(GetTonemappedBuffer(), BUFFER_PIXEL_FORMAT_DEFAULT, true);

  {
    static unsigned frame_count = 0;
    sutil::displayFps( frame_count++ );
  }
  
  frame_number++;

  glutSwapBuffers();
}

void glutKeyboardPress( unsigned char k, int x, int y )
{
    switch( k )
  {
    case '1':
    {
      show_buffer = 0;
      break;
    }
    case '2':
    {
      show_buffer = 1;
      break;
    }
    case '3':
    {
      show_buffer = 2;
      break;
    }
    case '4' :
    {
      show_buffer = 3;
      break;
    }
    case '5':
    {
      show_buffer = 4;
      break;
    }
    case( 'q' ):
    case( 27 ): // ESC
    {
      DestroyOptixContext();
      exit(0);
    }
  }

  glutPostRedisplay();
}

void glutMousePress( int button, int state, int x, int y )
{
  if( state == GLUT_DOWN )
  {
    mouse_button = button;
    mouse_prev_pos = optix::make_int2( x, y );
  }
  else
  {
    // nothing
  }

  glutPostRedisplay();
}

void glutMouseMotion( int x, int y)
{
  if( mouse_button == GLUT_RIGHT_BUTTON )
  {
    const float dx = static_cast<float>( x - mouse_prev_pos.x ) /
                     static_cast<float>(width );
    const float dy = static_cast<float>( y - mouse_prev_pos.y ) /
                     static_cast<float>(height );
    const float dmax = fabsf( dx ) > fabs( dy ) ? dx : dy;
    const float scale = fminf( dmax, 0.9f );
    camera_eye = camera_eye + (camera_lookat - camera_eye)*scale;
    camera_changed = true;
  }
  else if( mouse_button == GLUT_LEFT_BUTTON )
  {
    const optix::float2 from = { static_cast<float>(mouse_prev_pos.x),
                                 static_cast<float>(mouse_prev_pos.y) };
    const optix::float2 to   = { static_cast<float>(x),
                                 static_cast<float>(y) };

    const optix::float2 a = { from.x / width, from.y / height };
    const optix::float2 b = { to.x   / width, to.y   / height };

    camera_rotate = arcball.rotate( b, a );
    camera_changed = true;
  }

  mouse_prev_pos = optix::make_int2( x, y );

  glutPostRedisplay();
}

void glutResize( int w, int h )
{
  if ( w == (int)width && h == (int)height ) return;

  camera_changed = true;

  width = w;
  height = h;

  sutil::ensureMinimumSize(width, height);

  sutil::resizeBuffer(GetOutputBuffer(), width, height);
  sutil::resizeBuffer(GetTonemappedBuffer(), width, height);
  sutil::resizeBuffer(GetNormalBuffer(), width, height);
  sutil::resizeBuffer(GetAlbedoBuffer(), width, height);
  sutil::resizeBuffer(GetDenoisedBuffer(), width, height);

  glViewport(0, 0, w, h);
  postprocessing_must_init = true;

  glutPostRedisplay();
}

int main( int argc, char** argv )
{
  try
  {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_ALPHA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(width, height);
    glutInitWindowPosition(100, 100);
    glutCreateWindow(SAMPLE_NAME);
    glutHideWindow();

    glewInit();

    // load the ptx source associated with tutorial number
    std::string tutorial_ptx_path("raymarching.cu");
    tutorial_ptx = sutil::getPtxString( SAMPLE_NAME, tutorial_ptx_path.c_str() );

    CreateOptixContext();
    CreateTextureSamplers();
    CreateSceneGeometry();

    // Let's setup our camera
    camera_eye = optix::make_float3(512.0f, 512.0f, 512.0f);
    camera_lookat = optix::make_float3(0.0f, 0.0f, 0.0f);
    camera_up = optix::make_float3(0.0f, 1.0f, 0.0f);

    camera_rotate = optix::Matrix4x4::identity();
    
    context->validate();

    // Initialize GL state
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 1, 0, 1, -1, 1);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glViewport(0, 0, width, height);

    glutShowWindow();
    glutReshapeWindow(width, height);

    // register glut callbacks
    glutDisplayFunc(glutDisplay);
    glutIdleFunc(glutDisplay);
    glutReshapeFunc(glutResize);
    glutKeyboardFunc(glutKeyboardPress);
    glutMouseFunc(glutMousePress);
    glutMotionFunc(glutMouseMotion);
    glutCloseFunc(DestroyOptixContext);  // this function is freeglut-only

    glutMainLoop();

  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}