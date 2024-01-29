/*
  PhinGen3 - UGen for SuperCollider
  new version with static physical parameters

  A series of vertically moving masses connected with springs is used
  to generate periodically updated waveforms. 

  VERSION 3: The Waveform is updated according to a given modulatable frequency
 

  copyright 2011 by Luc Doebereiner (luc.doebereiner@gmail.com)
	
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include "SC_PlugIn.h"
#include <stdlib.h>

static InterfaceTable *ft;

struct mass {
  float posOld1;
  float posOld2;
  float Xinit;
  float force;
  float mass;
  float dX;
  float minX;
  float maxX;
  float posNew;
  float speedNew;
  float forceOut;  
};

struct lia {
  float stiffness;
  float viscosity; 
  float D2;
  float length;
  float distanceOld;
  float pos1;
  float pos2;
  float posOld1;
  float posOld2;
  float force1;
  float force2;
  float Lmin;
  float Lmax;
};


inline void mass_init(mass *m, float mass) {
  m->posOld1 = 0;
  m->posOld2 = 0;
  m->Xinit = 0;
  m->force = 0;
  //  m->mass = 50000;
  m->mass = mass;
  m->dX = 0;
  m->minX = -100000;
  m->maxX = 100000;
  m->posNew = 0;
  m->speedNew = 0;
  m->forceOut = 0;  
}

inline void mass_next(mass *m) {
  if (m->mass > 0)
    m->posNew = m->force/m->mass + 2*m->posOld1 - m->posOld2;
  else m->posNew = m->posOld1;
  m->posNew = std::max(std::min(m->maxX, m->posNew), m->minX);
  m->posNew += m->dX;
  m->posOld1 += m->dX; 
  
  m->speedNew = m->posOld1 - m->posOld2;

  m->posOld2 = m->posOld1;
  m->posOld1 = m->posNew;

  m->force = 0;
  m->dX = 0;
}

inline void mass_addForce(mass *m, float addForce) {
  m->force += addForce;
}

inline void mass_set_mass(mass *m, float mass) {
  m->mass = mass;
}

inline float mass_getPosition(mass *m) {
  return m->posNew;
}


inline void mass_setPosition(mass *m, float position) {
  // clearing history for stability
  m->posOld1 = position;
  m->posOld2 = position;
  m->force = 0;
  m->posNew = position;
}



inline void link_init(lia *l, float stiffness, float viscosity, float d2) {
  //  l->stiffness = 200;
  l->stiffness = stiffness;
  //  l->viscosity = 200;
  l->viscosity = viscosity;
  //  l->D2 = 20;
  l->D2 = d2;
  l->length = 0;
  l->distanceOld = 0;
  l->pos1 = 0;
  l->pos2 = 0;
  l->posOld1 = 0;
  l->posOld2 = 0;
  l->force1 = 0;
  l->force2 = 0;
  l->Lmin = 0;
  l->Lmax = 1000;
}


inline void link_setPos1(lia *l, float p1) {
  l->pos1 = p1;
}

inline void link_setPos2(lia *l,float p2) {
  l->pos2 = p2;
}

inline void link_set_stiffness(lia *l,float stiffness) {
  l->stiffness = stiffness;
}

inline void link_set_viscosity(lia *l,float viscosity) {
  l->viscosity = viscosity;
}

inline void link_set_d2(lia *l, float d2) {
  l->viscosity = d2;
}

inline float link_getForce1(lia *l) {
  return l->force1;
}

inline float link_getForce2(lia *l) {
  return l->force2;
}


inline void link_next(lia *l) {
  float distance;

  distance = (l->pos2 - l->pos1);
  if (distance<0) distance = -distance;
  l->force1 =  l->stiffness*(distance-(l->length)) + l->viscosity*(distance - l->distanceOld);
  l->distanceOld = distance;

  if (distance > l->Lmax) l->force1=0;
  if (distance < l->Lmin) l->force1=0;
  
  if (distance != 0) {
    l->force1 = l->force1 * (l->pos2 - l->pos1) / distance;
  }

  l->force2 = -l->force1 + (l->posOld2 - l->pos2)*l->D2;
  l->force1 += (l->posOld1 - l->pos1)*l->D2;

  l->posOld1 = l->pos1;
  l->posOld2 = l->pos2;
}

inline float fl_mod(float x, int m) {
  int truncated_x = x;
  float dec_part = x - truncated_x;
  return (truncated_x % m) + dec_part;
}


struct PhinGen3 : public Unit
{
  mass *masses; // the array of masses
  lia *links; // the array of links
  //  float *saved_triggers; // all the triggers, which have been received during the previous cycle are "saved" here
  float* phases; // the phases in fractional indices into the mass and link arrays
  float* sample_incs; // the sample increments for the readers
  int nreaders; // the number of "readers"/outputs
  int length; // the size of the mass and link arrays
  float mmass; // the mass of the masses
  float stiffness; // the stiffness of the links
  float viscosity; // the viscosity of the links
  float d2; // damping of the links
  //  bool spring_update_necessary; // if a period has been finished an update is required
  int distance;
  int update_duration; // duration (in smpls) between updates
  int update_counter; // if = update_duration, we need to update the springs
  float *last_positions; // last positions of the masses, used to interpolate
};



extern "C"
{
  void load(InterfaceTable *inTable);
  void PhinGen3_next_aa(PhinGen3 *unit, int inNumSamples);
  void PhinGen3_Ctor(PhinGen3* unit);	
  void PhinGen3_Dtor(PhinGen3* unit);	
  // insert all the necessary functions in here
};

void PhinGen3_Ctor(PhinGen3* unit) 
{
	
  // if (INRATE(0) == calc_FullRate) {
   //   if (INRATE(1) == calc_FullRate) {
  //       SETCALC(PhinGen3_next_aa);
  //   }
  //   else {
  //     SETCALC(PhinGen3_next_ak);
  //   }
  // }
  // else if (INRATE(1) == calc_FullRate) {
  //   SETCALC(PhinGen3_next_ka);
  // }
  // else {
  //    SETCALC(PhinGen3_next_kk);
  // }

  SETCALC(PhinGen3_next_aa);


  unit->length = sc_clip(IN0(5),2,300);
  unit->distance = sc_clip(IN0(6),1,300);
  //  printf("size: %i",unit->length * unit->distance); fflush(stdout);
  // * 100 becuase it is the maximum distance,  now 300 
  // * 50 because it is the maximum length, now 300
  unit->masses = (mass*)RTAlloc(unit->mWorld, (300 * 300 * sizeof(struct mass)));
  unit->links = (lia*)RTAlloc(unit->mWorld, (300 * 300 * sizeof(struct lia)));
  unit->last_positions = (float*)RTAlloc(unit->mWorld, (50 * 100 * sizeof(float)));
  //  unit->saved_triggers = (float*)RTAlloc(unit->mWorld, (unit->length * sizeof(float)));
  unit->nreaders = IN0(0);
  unit->phases = (float*)RTAlloc(unit->mWorld, (IN0(0) * sizeof(float)));
  unit->sample_incs = (float*)RTAlloc(unit->mWorld, (IN0(0) * sizeof(float)));
  unit->mmass = IN0(4); 
  unit->stiffness = IN0(1);
  unit->viscosity = IN0(2);
  unit->d2 = IN0(3);
  unit->update_duration = unit->mRate->mSampleRate / sc_clip(IN0(7),0.5,15000); 
  unit->update_counter = 0;

  // init the masses and links
  for (int i = 0; i < (unit->length * 300); i++) {
    mass_init(&(unit->masses[i]), unit->mmass);
    link_init(&(unit->links[i]), unit->stiffness, unit->viscosity, unit->d2);
  }

  for (int i = 0; i < unit->nreaders; i++) {
    unit->phases[i] = 0.0;
  }

  //  printf("did init"); fflush(stdout);

  // for (int i = 0; i < unit->length; i++) {
  //   unit->saved_triggers[i] = 0.0;
  // }

  // unit->tg_smpl_inc = unit->tg_bp_vals[(unit->tg_max_idx)-1] / 
  //   ((float) unit->mRate->mSampleRate / sc_max(unit->mInBuf[0][0],.001f)); 
  PhinGen3_next_aa(unit, 1);
}


void PhinGen3_Dtor(PhinGen3 *unit)
{
	RTFree(unit->mWorld, unit->masses);
	RTFree(unit->mWorld, unit->links);
	RTFree(unit->mWorld, unit->last_positions);
	RTFree(unit->mWorld, unit->phases);
	RTFree(unit->mWorld, unit->sample_incs);
	//	RTFree(unit->mWorld, unit->saved_triggers);
}

inline float lin_interpolation(float x, float y1, float y2) 
{
  return y1 + ((y2 - y1) * x);
}


void PhinGen3_next_aa(PhinGen3 *unit, int inNumSamples)
{  
  //  printf("calling next\n"); fflush(stdout);
  // copy values from the unit

  float stiffness = unit->stiffness;
  float mmass = unit->mmass;
  float viscosity = unit->viscosity;
  float d2 = unit->d2;
  int length = unit->length;
  int distance = unit->distance;
  int nmasses = length * distance;
  mass *masses = unit->masses;
  lia *links = unit->links;
  float* phases = unit->phases;
  int nreaders = unit->nreaders;
  //  float *saved_triggers = unit->saved_triggers;
  //  bool spring_update_necessary = unit->spring_update_necessary;
  int update_counter = unit->update_counter;
  int update_duration = unit->update_duration;
  float update_x; //counter / duration;
  float *last_positions = unit->last_positions;
  float *sample_incs = unit->sample_incs;


  int x0, x1, x2, x3; // the indices for the cubic interpolation


  // input and output buffers
  //  float *out = OUT(0);
  //  float *freq_in0 = IN(0); // freq
  //  float *stiffness_in1 = IN(1); // stiffness
  //  float *viscosity_in2 = IN(2); // viscosity
  //  float *d2_in3 = IN(3); // d2
  //  float *mass_in4 = IN(4); // mass
  float *length_in5 = IN(5); // number of controllable breakpoints
  float *distance_in6 = IN(6); // distances between them
  float *update_freq = IN(7);
  float **in = unit->mInBuf; // the array of frequencies and input triggers


	
  // the audio loop
  for (int i=0; i < inNumSamples; ++i) {
    //    printf("looping\n"); fflush(stdout);
    length = sc_clip(length_in5[i],2,300);
    distance = sc_clip(distance_in6[i],1,300);
    nmasses = length * distance; // total number of masses

    //    printf("setting sample incs\n"); fflush(stdout);
    for (int k = 0; k < nreaders;k++) {
      sample_incs[k] = nmasses * in[8+k][i] / float(unit->mRate->mSampleRate);
    }
    //    printf("sample incs set\n"); fflush(stdout);

    if (update_counter >= update_duration) {
      //      printf("updating\n"); fflush(stdout);
      // update link and mass parameters here

      // copy current values to last positions before updating 
      for (int last_pos_idx = 0; last_pos_idx < nmasses; last_pos_idx++) {
	last_positions[last_pos_idx] = mass_getPosition(&masses[last_pos_idx]);
      }

      //      mmass = mass_in4[i];
      //      stiffness = stiffness_in1[i];
      //      viscosity = viscosity_in2[i];
      //      d2 = d2_in3[i];
      

      // every nth mass acts as a controllable breakpoint
      for (int idx_bp = 0; idx_bp < length; idx_bp++) {
	
	mass_setPosition(&masses[idx_bp*distance], in[8+nreaders+idx_bp][i]);
	//	printf("breakpoint:%f\n",in[6+idx_bp][i]); fflush(stdout);
      }


      // updating the springs
      for (int idx = 0; idx < nmasses; idx++) {

	// setting the parameters
	mass_set_mass(&masses[idx], mmass);
	//	printf("setting mass to:%f\n",mmass);fflush(stdout);
	link_set_stiffness(&links[idx], stiffness);
	link_set_viscosity(&links[idx], viscosity);
	link_set_d2(&links[idx], d2);

	// calculating the link before
	if (idx >= 1) 
	  link_setPos1(&links[idx], mass_getPosition(&masses[idx-1]));
	else
	  link_setPos1(&links[idx], mass_getPosition(&masses[nmasses-1]));
	link_setPos2(&links[idx], mass_getPosition(&masses[idx]));
	link_next(&links[idx]);
	// setting the additional force for the previous mass
	if (idx >= 1)
	  mass_addForce(&masses[idx-1], link_getForce1(&links[idx]));
	else
	  mass_addForce(&masses[nmasses-1], link_getForce1(&links[idx]));
	// setting the additional force for this mass
	mass_addForce(&masses[idx], link_getForce2(&links[idx]));
	mass_next(&masses[idx]); // calculating this mass
      };

      update_duration = unit->mRate->mSampleRate / sc_clip(update_freq[i],0.5,15000);
      update_counter = 0;

    }

    update_x = (float) update_counter / (float) update_duration;
    update_counter += 1;

    //    printf("setting phases and output, nreaders: %i\n", nreaders); fflush(stdout);
 
    //    out0[i] = 0.0;
    //out1[i] = 1.0;

   for (int k = 0;k < nreaders;k++) {
      //      printf("setting phases %i\n", k); fflush(stdout);
      phases[k] = fl_mod(phases[k] + sample_incs[k], nmasses); // updating the phase
      // getting the four surrounding indices for the cubic interpolation
      x1 = int(phases[k]);
      x0 = (x1 == 0) ? (nmasses - 1) : (x1 - 1); // wrap around
      x2 = (x1 + 1) % nmasses;
      x3 = (x2 + 1) % nmasses;
      (unit->mOutBuf[k])[i] = cubicinterp(phases[k] - x1, 
    			 lin_interpolation(update_x, last_positions[x0], mass_getPosition(&masses[x0])),
    			 lin_interpolation(update_x, last_positions[x1], mass_getPosition(&masses[x1])),
    			 lin_interpolation(update_x, last_positions[x2], mass_getPosition(&masses[x2])),
    			 lin_interpolation(update_x, last_positions[x3], mass_getPosition(&masses[x3])));
    }
    // //    printf("done with output\n"); fflush(stdout);

  }
  
  // store values in the unit
  unit->stiffness = stiffness;
  unit->mmass = mmass;
  unit->viscosity = viscosity;
  unit->d2 = d2;
  unit->length = length;
  unit->masses = masses;
  unit->links = links;
  unit->phases = phases;
  //  unit->saved_triggers = saved_triggers;
  //  unit->spring_update_necessary = spring_update_necessary;
  unit->distance = distance;
  unit->update_counter = update_counter;
  unit->update_duration = update_duration;
  unit->last_positions = last_positions;
}


PluginLoad(PhinGen3)
{
  ft = inTable;

  DefineDtorUnit(PhinGen3);
}


// void load(InterfaceTable *inTable)
// {
//   ft = inTable;

//   DefineDtorUnit(PhinGen3);
// }


