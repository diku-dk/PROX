#version 410 core

struct LightInfo
{
  vec3  position; // Position of light in world
  vec3  target;   // Position of light target in world
  float cutoff_angle;
  float attenuation;
  vec3  Is;       // specular color
  vec3  Id;       // diffuse color
  vec3  Ia;       // ambient color
};

struct MaterialInfo
{
  vec3 Ks;                 // specular color
  vec3 Kd;                 // diffuse color
  vec3 Ka;                 // ambient color
  float specular_exponent; // specular 'power'
  bool use_texture;        // True if Kd should be replaced with solid texture
  float alpha;             // transparent (0.0f) to opaque (1.0f)
};

in vec3 position_eye;        // position in camera eye space
in vec3 normal_eye;          // normal in camera eye space
in vec3 tex_coord;           // texture coordinate
in vec4 position_light[4];   // position in light coordinate frame

out vec4 frag_color;

uniform mat4              view_matrix;
uniform sampler3D         solid_texture;
uniform MaterialInfo      material;
uniform int               number_of_lights;
uniform LightInfo         lights[4];

/**
 * @param P  surface point position in Eye frame
 * @param N  surface point normal in Eye frame
 */
vec3 compute_light_intensity(
                            int light_source
                           , MaterialInfo model
                           , vec3 P
                           , vec3 N
                           )
{
  LightInfo light = lights[light_source];

  // Light position in eye frame
  vec3 S = vec3 (view_matrix * vec4 (light.position, 1.0));

  // light direction in eye frame
  vec3 D = vec3 (view_matrix * vec4 (normalize(light.target - light.position), 0.0));

  float distance = length(S - P);

  vec3 L = normalize (S - P);

  //--- Ambient intensity ------------------------------------------------------
  vec3 La = light.Ia * model.Ka;

  //--- Diffuse intensity ------------------------------------------------------
  float diffuse_factor =  max (min ( dot (L, N), 1.0f), 0.0f);

  vec3 Ld = light.Id * model.Kd * diffuse_factor;

  //---- Specular intensity ----------------------------------------------------
  vec3 R = reflect (-L, N);
  vec3 V = normalize (-P);

  float specular_factor = pow ( max( dot (R, V), 0.0) , model.specular_exponent);

  vec3 Ls = light.Is * model.Ks * specular_factor;

  //--- Splot light cone and shadows

  float attenuation = 1.0 / (1.0 + light.attenuation * pow(distance, 2));

  //--- Final intensity ------------------------------------------------------------
  return ( attenuation*(Ls + Ld) + La );
}

void main()
{
  MaterialInfo model = material;

  vec3 Kd_texture =  texture( solid_texture, tex_coord ).rgb;
  model.Kd        = material.use_texture ? Kd_texture : material.Kd;

  vec3 accum = vec3(0.0,0.0,0.0);

  for (int i=0; i<number_of_lights; ++i)
  {
    vec3 color = compute_light_intensity(
                                    i
                                  , model
                                  , position_eye
                                  , normal_eye
                                  );
    accum = accum + color;
  }

  frag_color = vec4 (accum, model.alpha);
}

