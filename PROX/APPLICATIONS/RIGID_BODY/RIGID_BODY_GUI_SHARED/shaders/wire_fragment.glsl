#version 150

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

in vec3 v_position_eye;
in vec3 v_normal_eye;
in vec3 v_tex_coord;
in vec3 v_bary_coord;

out vec4 frag_color;

uniform vec3         wire_color;
uniform mat4         view_matrix;
uniform sampler3D    solid_texture;
uniform int          number_of_lights;
uniform LightInfo    lights[4];
uniform MaterialInfo material;

/**
 * @param P  surface point position in Eye frame
 * @param N  surface point normal in Eye frame
 */
vec3 compute_light_intensity(LightInfo light, MaterialInfo model, vec3 P, vec3 N)
{
  // Light position in eye frame
  vec3 S = vec3 (view_matrix * vec4 (light.position, 1.0));

  // light direction in eye frame
  vec3 D = vec3 (view_matrix * vec4 (normalize(light.target - light.position), 0.0));

  float distance = length(S - P);

  vec3 L = normalize (S - P);

  float attenuation = 1.0 / (1.0 + light.attenuation * pow(distance, 2));

  float light_angle = degrees( acos(dot(-L, D)));

  if(light_angle > light.cutoff_angle )
  {
    attenuation = 0.0;
  }

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

  //--- Final intensity ------------------------------------------------------------
  return ( attenuation*(Ls + Ld) + La );
}


float wire_edge_factor()
{
  vec3 d = fwidth(v_bary_coord);

  float thickness = 2.0; // 2014-10-06 Kenny: This could be a uniform providing greather control?

  vec3 a3 = smoothstep(vec3(0.0), d*thickness, v_bary_coord);

  return min(min(a3.x, a3.y), a3.z);
}

void main()
{
  MaterialInfo model = material;

  vec3 Kd_texture =  texture( solid_texture, v_tex_coord ).rgb;

  model.Kd        = material.use_texture ? Kd_texture : material.Kd;

  vec3 accum = vec3(0.0,0.0,0.0);
  for (int i=0; i<number_of_lights; ++i)
  {
    vec3 color = compute_light_intensity(lights[i], model, v_position_eye, v_normal_eye);
    accum = accum + color;
  }

  frag_color = vec4( mix(wire_color, accum, wire_edge_factor()), 1.0 );
//frag_color = vec4( vec3(edge_factor()), 1.0 );
}

