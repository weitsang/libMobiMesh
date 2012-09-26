/*
 *  StatusInfo.h
 *  MobiMesh
 *
 *  Created by Ong Yuh Shin on 3/10/11.
 *  Copyright 2011 National University of Singapore. All rights reserved.
 *
 */

#ifndef __STATUS_INFO_H_
#define __STATUS_INFO_H_

namespace MobiMesh {

enum StatusBits {
    DELETED               = 1,    // Item has been deleted
    LOCKED                = 2,    // Item is locked.
    SELECTED              = 4,    // Item is selected.
    HIDDEN                = 8,    // Item is hidden.
    FEATURE               = 16,   // Item is a feature or belongs to a feature.
    TAGGED                = 32,   // Item is tagged.
    TAGGED2               = 64,   // Alternate bit for tagging an item.
    FIXEDNONMANIFOLD      = 128,  // Item was non-two-manifold and had to be fixed
    UNUSED                = 256   // Unused
};

class StatusInfo {
	
public:
    typedef unsigned int value_type; 
    StatusInfo() : status_(0) {}
	
    bool     deleted()  const           { return is_bit_set(DELETED); }
    void set_deleted (bool _b)          { change_bit(DELETED, _b); }
    bool     locked()   const           { return is_bit_set(LOCKED); }
    void set_locked  (bool _b)          { change_bit(LOCKED, _b); }
    bool     selected() const           { return is_bit_set(SELECTED); }
    void set_selected(bool _b)          { change_bit(SELECTED, _b); }
    bool     hidden()   const           { return is_bit_set(HIDDEN); }
    void set_hidden  (bool _b)          { change_bit(HIDDEN, _b); }
    bool     feature()  const           { return is_bit_set(FEATURE); }
    void set_feature (bool _b)          { change_bit(FEATURE, _b); }
    bool     tagged()  const            { return is_bit_set(TAGGED); }
    void set_tagged  (bool _b)          { change_bit(TAGGED, _b); }
    bool     tagged2() const            { return is_bit_set(TAGGED2); }
    void set_tagged2 (bool _b)          { change_bit(TAGGED2, _b); }
    bool     fixed_nonmanifold() const  { return is_bit_set(FIXEDNONMANIFOLD); }
    void set_fixed_nonmanifold(bool _b) { change_bit(FIXEDNONMANIFOLD, _b); }
	
    unsigned int bits() const             { return status_; }
    void     set_bits(unsigned int _bits) { status_ = _bits; }
	
    bool is_bit_set(unsigned int _s) const    { return (status_ & _s) > 0; }
    void    set_bit(unsigned int _s)          { status_ |= _s; }
    void  unset_bit(unsigned int _s)          { status_ &= ~_s; }
    void change_bit(unsigned int _s, bool _b) { if (_b) status_ |= _s; else status_ &= ~_s; }
	
private:
	value_type status_;
};

}

#endif
