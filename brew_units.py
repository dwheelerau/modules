# -*- coding: iso-8859-1 -*-
# Homebrewing recipe sharing and calculation software
# copyright (c) 2004 paul sorenson
# $Id: units.py 400 2006-03-10 23:25:02Z  $

# This file is part of Brewsta
#
# Brewsta is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Brewsta is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Brewsta if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import logging
# NOTE about logging.  Unit objects are called in the grid refresh which gets
# called many times and is performance sensitive.  That is why logging is
# commented out for now.  It will work with logging turned on, just noticeably
# slower refreshes.
log = logging.getLogger('brewsta')
##glog = logging.getLogger('grid')

import re

united = re.compile('(-?\d+\.?\d*) *(.*?)$')

DEF_FORMAT = "%(val)2.1f %(unit)s"

USGAL_PER_LITRE = 0.26417205
USLB_PER_KG = 2.2046226
IMPLB_PER_KG = 1.8357323

import weakref
allUnits = weakref.WeakKeyDictionary()
def toSi():
    for u in allUnits.keys():
        u.toSi()
def fromSi():
    for u in allUnits.keys():
        u.fromSi()

# Maps references to units
unitsRefs = weakref.WeakValueDictionary()
def getUnitByRef(ref):
    try:
        return unitsRefs[ref]
    except:
        return None

def unitDefaults():
    urefs = unitsRefs.keys()
    urefs.sort()
    for uref in urefs:
        yield uref, unitsRefs[uref].getUnits()

# Unit systems
UK = 'UK'
US = 'US'
SI = 'SI'

# When switching default units (cf systems), these are the units that are
# most likely to correspond.
FROM_SI = {'g': 'oz', 'kg': 'lb', 'L': 'gal', 'C': 'F', 'HWE': 'PPG', 'EBC': 'SRM', 'L/kg': 'qt/lb', 'PL': 'PG'}
# Add reverse mapping
TO_SI = {}
for k in FROM_SI.keys(): TO_SI[FROM_SI[k]] = k

class NonLinearUnitDesc(object):
    '''
    Generalized unit conversion

    Using the get/set has the side effect of changing the owners current
    units.
    '''
    def __init__(self, toBase, fromBase, foreignUnit = '', format = None):
        '''
        Constructor

        toBase: used when setting native units from foreign units.
        fromBase: used when getting foriegn units from native units.
        foreignUnit: name of the unit.  Should be the same as the name of the
        attribute you are creating when lowercased.
        format: format string for displaying foreign unit.
        '''
        self.toBase = toBase
        self.fromBase = fromBase
        self.foreignUnit = foreignUnit
        self.format = format

    def __get__(self, obj, type=None):
        '''
        Get the unit converted to a foreign unit.
        '''
        if obj.value != None:
            foreignUnit = self.fromBase(obj.value)
        else:
            foreignUnit =  None
        if self.format:
            format = self.format
        else:
            format = obj.format
        return foreignUnit, self.foreignUnit, format

    def __set__(self, obj, value):
        if value != None:
            obj.value = self.toBase(value)


class LinearUnitDesc(NonLinearUnitDesc):
    '''
    Implement descriptors for linear unit conversion.
    Useful for creating attributes in sub classes of Unit() eg:
        g = LinearUnitDesc(0.001, 'g') where the class represents mass in kg.
    '''
    def __init__(self, factor = 1.0, foreignUnit = '', format = None):
        '''
        factor: multiplier used when setting native units from foreign units.
        foreignUnit: name of the unit.  This should be the same as the name
        of the attribute you are creating (case insensitive).
        format: format string for displaying foreign unit.
        '''
        NonLinearUnitDesc.__init__(self,
             lambda f: float(factor) * f,
             lambda n: n / float(factor),
             foreignUnit,
             format
             )


class SystemAlias(object):
    '''
    A way to alias units which might change depending on preferences.
    Eg pounds = SystemAlias({UK: brpounds, US: uspounds})

    Changing the classes "system" attribute is a method for making wholesale
    unit changes at runtime (at least SI Imperial and US).  Note that this has
    no effect on "default" units, it just changes the definition of say "pound"
    between imperial and US.
    '''
    def __init__(self, map):
        self.map = map

    def __get__(self, obj, type=None):
        return self.map[obj.system].__get__(obj)

    def __set__(self, obj, value):
        self.map[obj.system].__set__(obj, value)


class Alias(object):
    def __init__(self, existingUnit):
        self.unit = existingUnit

    def __get__(self, obj, type=None):
        return self.unit.__get__(obj)

    def __set__(self, obj, value):
        self.unit.__set__(obj, value)


class Unit(object):
    '''
    Value class for a unit.
    Self.value always contains base (native) units.  This is not the same
    as default units.
    WARNING: attributes should be defined as lower case strings if setFromString
    is used since it lowercases all strings.  TODO: review this.  Check
    setFromString() for latest.
    '''
    def __init__(self, nativeUnit = '', value = None, format = DEF_FORMAT, aliases = None, ref = None):
        '''
        nativeUnit: text string to represent native unit, eg 'kg'.
        value: initial value in native units.
        format: default format used when returning string representation.  The format
        string will get interpolated from a map with keys: 'val', 'unit', 'format'.
        aliases: provides a mechanism to lookup arbitrary units strings that would
        otherwise not be legal as attributes, eg '%'.
        '''
        self.currentUnit = self.nativeUnit = nativeUnit
        self.value = value
        self.set(value)
        self.format = format
        self.aliases = aliases
        allUnits[self] = None
        self.ref = ref
        if ref:
            unitsRefs[ref] = self

    def setAliases(self, aliases):
        self.aliases = aliases

    def unAlias(self, unit):
        '''
        Looks up possible alias for a particular unit.
        If not found, return the orginal unit
        '''
        try:
            result = self.aliases[unit]
        except:
            result = unit
        return result

    def getUnits(self):
        return self.currentUnit

    def setUnits(self, units = None):
        '''
        Set the current units.  Values will be retrieved in these units
        unless specified otherwise.

        units: specifies current units for certain access.  If units is not
        specified then sets units to native units.
        '''
        if not units:
            units = self.nativeUnits
        self.currentUnit = units

    def get(self, units = None):
        '''
        Get current value as (value, units).

        Value is returned in the units specified or the current units otherwise.
        '''
        if not units:
            units = self.currentUnit
        units = self.unAlias(units)
##        log.debug('Using units: "%s"' % units)
        return getattr(self, units)

    def getNativeValue(self):
        return self.value

    def setNativeValue(self, value):
        self.value = value

    def set(self, value, units = None):
        '''
        Set new value, optionally specifying what units the value is in.

        units: if set, the value is converted assuming these units.  If it
        is None then current units are assumed.
        '''
        if not units:
            units = self.currentUnit
        units = self.unAlias(units)
##        log.debug('setting units: (value %s)(units %s)' % (value, units))
        # You would think it makes sense to do a hasattr test right here to
        # avoid writing garbage attributes to the object.  However it seems
        # that hasattr fails for the first time.
        setattr(self, units, value)

    def setFromString(self, s):
        '''
        Parse a string with number possibly followed by units spec.

        For convenience, it returns the native value.  This could be used
        for parsing data files.
        '''
        if not s:
            return None
        m = united.match(s)
        def twoGroups(val, u = None):
            return float(val), u
        val, u = twoGroups(*m.groups())
##        glog.debug('setFromString: (text %s)(val %s)(units %s)' % (s, val, u))
        self.set(val, u)
        return self.getNativeValue()

    def __call__(self, s):
        return self.setFromString(s)

    def setSystem(klass, system):
        '''
        Set the measurement system globally for all instances.
        '''
        klass.system = system
    setSystem = classmethod(setSystem)

    def switch(self, map):
        '''
        If there is a mapping found for the currentUnit in map then
        switch the currentUnit to it.
        '''
        try:
            self.currentUnit = map[self.currentUnit]
        except:
            pass

    def fromSi(self):
        self.switch(FROM_SI)

    def toSi(self):
        self.switch(TO_SI)
##        self.currentUnit = self.nativeUnit

    def __str__(self):
        '''
        Get the string representation of the Unit in currentUnit.

        Consider taking a closer look at this with the profiler, this gets
        called many times when grids are refreshed.
        '''
        t = self.get()
        if t == None:
            return ''
##        log.debug("formatting units: '%s'" % str(t))
        map = dict(zip(('val', 'unit', 'format'), t))
##        log.debug(str(map))
        try:
            s = map['format'] % map
##            log.debug(s)
        except:
            s = map['unit']
        return s

    def unitGen(self):
        '''
        Generate list of descriptors - might be useful in enumerating
        units provided by an instance.
        '''
        d = self.__class__.__dict__
        for key in d.keys():
            if isinstance(d[key], NonLinearUnitDesc):
                yield key

    def sortedUnits(self):
        units = [unit for unit in self.unitGen()]
        units.sort()
        return units

Unit.setSystem(US)


class Unity(Unit):
    '''
    A convenience Unit subclass for dimensionless or single unit values.

    It basically leverages the formatting of the Unit class.
    '''

    def __init__(self, value = None, format = '%(val)2.2f', ref = None):
        Unit.__init__(self, 'fraction', value, format, ref=ref)

    # Use default format so it is convenient to override
    fraction = LinearUnitDesc(1, 'fraction')


class Percent(Unit):
    '''
    Format numbers as fractions or percentage.  Could subclass from Unity
    but would have to do some more trickery to show super class attributes
    in menu.
    '''
    def __init__(self, value = None, format = DEF_FORMAT, ref = None):
        Unit.__init__(self, '%', value, format, aliases = {'%': 'percent'}, ref=ref)

    fraction = LinearUnitDesc(1, 'fraction', format = '%(val)2.2f')
    percent = LinearUnitDesc(0.01, 'percent', format = '%(val)2.1f%%')


class Mass(Unit):
    def __init__(self, value = None, format = DEF_FORMAT, ref = None):
        Unit.__init__(self, 'kg', value, format, aliases = {'#': 'lb'}, ref=ref)

    kg = LinearUnitDesc(1.0, 'kg', format = "%(val)2.2f %(unit)s")
    tonne = LinearUnitDesc(1000, 't', format = "%(val)2.3f %(unit)s")
    g = LinearUnitDesc(0.001, 'g')
    gram = Alias(g)

    # US units
    lb = LinearUnitDesc(0.45359237, 'lb', format = "%(val)2.2f %(unit)s")
    pound = Alias(lb)
    oz = LinearUnitDesc(0.028349523, 'oz', format = "%(val)2.2f %(unit)s")
    ounce = Alias(oz)

    imppound = LinearUnitDesc(0.54474172, 'lb')
    impounce = LinearUnitDesc(0.034046358, 'oz', format = "%(val)2.2f %(unit)s")

##    lb = pound = SystemAlias({UK: imppound, US: uspound})
##    oz = ounce = SystemAlias({UK: impounce, US: usounce})


class Volume(Unit):
    def __init__(self, value = None, format = DEF_FORMAT, ref = None):
        Unit.__init__(self, 'L', value, format, ref=ref)

    litre = LinearUnitDesc(1.0, 'L')
    L = Alias(litre)
    ml = LinearUnitDesc(0.001, 'mL')
    # US units
    gal = LinearUnitDesc(3.7854118, 'gal', format = "%(val)2.2f %(unit)s")
    gallon = Alias(gal)
    qt = LinearUnitDesc(0.94635295, 'qt', format = "%(val)2.2f %(unit)s")
    quart = Alias(qt)
    
    bbl = LinearUnitDesc(158.98729, 'bbl', format = "%(val)2.2f %(unit)s")
    barrel = Alias(bbl)

    impgal = LinearUnitDesc(4.54609, 'gal')
    impquart = LinearUnitDesc(1.1365225, 'qt')
##    gal = gallon = SystemAlias({UK: impgal, US: usgal})
##    qt = quart = SystemAlias({UK: impquart, US: usquart})


class Temperature(Unit):
    def __init__(self, value = None, format = DEF_FORMAT, ref = None):
        Unit.__init__(self, 'C', value, format, ref=ref)

    C = LinearUnitDesc(1.0, 'C')
    degC = Alias(C)
    F = NonLinearUnitDesc(lambda f: (f - 32) * 5.0/9.0, lambda c: c*9.0/5.0 + 32, 'F', format = "%(val)2.0f %(unit)s")
    degF = Alias(F)
    deg = SystemAlias({SI: C, UK: C, US: F})


from SG import sgToBrix1, brixToSg1
from ingredients import SUCROSE_PPG, SUCROSE_HWE, EXTRACT_PPG, EXTRACT_HWE, PPG2HWE, USLB_PER_KG

class Gravity(Unit):
    '''
    Defines the sugar content of wort.
    '''
    def __init__(self, value = None, format = DEF_FORMAT, ref = None):
        Unit.__init__(self, 'SG', value, format, aliases = {'%': 'pot'}, ref=ref)

    SG = LinearUnitDesc(1.0, 'SG', format = '%(val)4.3f %(unit)s')
    GU = NonLinearUnitDesc(lambda gu: gu / 1000.0 + 1,
        lambda sg: (sg - 1) * 1000, 'GU', format = '%(val)2.0f %(unit)s')
    pot = NonLinearUnitDesc(lambda pot: pot * SUCROSE_PPG / 100000.0 + 1,
        lambda sg: (sg - 1) * 100000.0/SUCROSE_PPG, 'pot', format = '%(val)2.1f %(unit)s')
    Plato = Brix = NonLinearUnitDesc(brixToSg1, sgToBrix1, 'Brix', format = '%(val)4.2f %(unit)s')
    deg = B = Alias(Brix)


class TotalGravity(Unit):
    '''
    Total gravity is points * volume
    '''
    def __init__(self, value = None, format = '%(val)4.0f %(unit)s', ref = None):
        Unit.__init__(self, 'kg', value, format, ref=ref)

    kg = LinearUnitDesc(1.0, 'kg', format = '%(val)2.2f %(unit)s')
    PL = LinearUnitDesc(1.0 / SUCROSE_HWE, 'PL')
    PG = LinearUnitDesc(3.7854118 / SUCROSE_HWE, 'PG')
    lb = LinearUnitDesc(1.0/USLB_PER_KG, 'lb', format = '%(val)2.2f %(unit)s')


class Potential(Unit):
    '''
    Defines the malt contribution per mass of grain to each unit of
    wort.
    '''
    def __init__(self, value = None, format = DEF_FORMAT, ref = None):
        Unit.__init__(self, 'HWE', value, format, aliases = {'%': 'potential'}, ref=ref)

    HWE = LinearUnitDesc(1.0, 'HWE', format = '%(val)3.0f %(unit)s')
    PPG = LinearUnitDesc(PPG2HWE, 'PPG', format = '%(val)2.0f %(unit)s')
    # TODO: potential is not always given as Course Grind, As Is but this is
    # the most common spec.
    potential = CGAI = LinearUnitDesc(SUCROSE_HWE/100, '%', format = '%(val)2.0f %(unit)s')


class Colour(Unit):
    def __init__(self, value = None, format = DEF_FORMAT, ref = None):
        Unit.__init__(self, 'EBC', value, format, aliases = {'°': 'deg', '°L': 'L'}, ref=ref)

    EBC = LinearUnitDesc(1.0, 'EBC', format = '%(val)2.0f %(unit)s')
    SRM = LinearUnitDesc(1.97, 'SRM', format = '%(val)2.0f %(unit)s')
    L = Lovibond = Alias(SRM)
    deg = SystemAlias({SI: EBC, UK: EBC, US: SRM})


class WaterToGrain(Unit):
    '''
    Handles mash thickness and loss to grain.

    TODO: allow combinations of units.
    '''
    def __init__(self, value = None, format = DEF_FORMAT, ref = None):
        Unit.__init__(self, 'L/kg', value, format, aliases = {'L/kg': 'Lkg', 'qt/lb': 'qtlb'}, ref=ref)
    # Underscores get dicked around when auto creating the menu to change
    # units dynamically.
    Lkg = LinearUnitDesc(1.0, 'L/kg', format = "%(val)2.2f %(unit)s")
    qtlb = LinearUnitDesc(2.0863511, 'qt/lb', format = "%(val)2.2f %(unit)s")


if __name__ == "__main__":
    m = Mass()
    print m
    m.set(1)
    print; print m

    mass = Mass(1)
    print mass, mass.get(), mass.get('kg')
    print mass.tonne, mass.get('tonne')
    print mass.uspound, mass.pound
    print mass.imppound, mass
    print mass.g, mass.gram
    print
    mass.uspound = 1
    print mass.currentUnit
    print mass
    print
    mass.lb = 2
    print mass
    print
    mass.set(1, 'kg')
    print mass.kg
    print
    mass.set(2)
    print mass.kg
    print mass.pound
    print
    mass2 = Mass(1.0)
    mass.setSystem(UK)
    print mass.pound
    print; mass.setFromString("3.0 kg"); print mass.kg
    print; mass.setFromString("4.5"); print mass.kg
    print; mass.setFromString("1.5 lb"); print mass

    print
    t = Temperature(100)
    print t, t.C
    print t.F
    t.F = 32
    print t, t.C
    print t.F

    print
    g = Gravity(1.040)
    print g.SG, g.GU, g.Brix
    print g
    for s in g.sortedUnits():
        print s

    print
    v = Volume()
    v.L = 100.0
    print v.L, v.gal
    v.setFromString("50 L")
    print v.L, str(v)
    v.fromSi()
    print v.L, str(v)

    print
    print FROM_SI, TO_SI

    print
    print allUnits.keys()

    print
    fromSi()
    print m, mass, mass2, t, g, v

    pot = Potential()
    pot.PPG = SUCROSE_PPG
    print pot.PPG
    print pot
