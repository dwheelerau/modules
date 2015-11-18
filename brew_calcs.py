# Homebrewing recipe sharing and calculation software
# copyright (c) 2004 paul sorenson
# $Id: brewcalcs.py 20 2006-07-09 12:13:44Z pms $

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

from math import *
import logging, logging.config
#log = logging.getLogger('brew')

CFG = 'brewsta.cfg'
log = logging.getLogger('brewsta')
logging.config.fileConfig(CFG)

PPG2HWE = 8.345

def utilizationTinseth(boilGravity, time):
    '''
    Utilization is a dimensionless value that expresses what fraction
    of the hops bittering capacity is used.
    '''
    fG = 1.65 * pow(0.000125, (boilGravity - 1))
    fT = (1 - exp(-0.04 * time)) / 4.15
    return fG * fT

def gravity(self, PPG, recipeVolume):
    return PPG() * PPG2HWE / self.recipeVolume

class Product:
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return self.name

    def tinsethByAA(self, boilGravity, time):
        '''
        This is a dimensionless characteristic of the product (ie hops) when
        the boil conditions are factored in.
        The factors that ultimately influence the beer bitterness include:
            Alpha acids in hops
            Amount of hops
            the timing of the hop addition, longer boil means more bitterness.
            final recipe volume
            boil gravity (affects utilization)

        IBU units are mg/litre
        '''
        return 0.0

    def PPG(self):  # points/lb/gal
        return 0.0

    def HWE(self):  # points/kg/litre
        return self.PPG() * PPG2HWE

class Grain(Product):
    def __init__(self, name, potential, colour):
        Product.__init__(self, name)
        self.potential = potential
        self.colour = colour

    def PPG(self):
        return self.potential

class Hops(Product):
    def __init__(self, name, AA):
        '''
        Specialized product that yields alpha acids for bittering.
        AA: alpha acids expressed as a fraction (%/100)
        '''
        Product.__init__(self, name)
        self.AA = AA    # AA is a fraction (not percentage value)

    def tinsethByAA(self, boilGravity, time):
        # Each gram of hops will bitter one litre this much
        return self.AA * utilizationTinseth(boilGravity, time)

class Ingredient:
    def __init__(self, product, amount, time = None):
        self.product = product
        self.amount = amount
        self.time = time

    def PG(self):
        # *** wrong *** amount is SI
        return self.amount * self.product.PPG()

##    def HWE(self):
##        return self.amount * self.product.PPG() * PPG2HWE

    def tinsethAcids(self, boilGravity):
        '''Total bittering in mg'''
        return self.amount * 1000000 * self.product.tinsethByAA(boilGravity, self.time)

    def __str__(self):
        return str('ingredient: "%s"  amount: %f  time: %s"' % (self.product, self.amount, self.time))

class Recipe:
    def __init__(self, recipeVolume = 23, ingredients = None):
        self.acids = self.tinsethAcids
        self.recipeVolume = recipeVolume
        self.boilVolume = recipeVolume
        self.ingredients = ingredients
        if not self.ingredients:
            self.ingredients = []
#        self.boilGravity = 1.000

    def addIngredient(self, ingredient):
        self.ingredients.append(ingredient)
        log.debug('added ingredient "%s"' % ingredient)

    def IBU(self):
        '''
        Calculate the IBU from contribute ingredient acids divided
        by recipe volume.  mg/litre
        '''
        acids = 0.0
        boilGravity = self.gravity()
        for ingredient in self.ingredients:
            acids += self.acids(ingredient, boilGravity)
        return acids / self.recipeVolume

    def tinsethAcids(self, ingredient, boilGravity):
        return ingredient.tinsethAcids(boilGravity)

    def PG(self):
        ''' Points per gallon (ie PPG * pounds) '''
        totalPg = 0.0
        for ingredient in self.ingredients:
            pg = ingredient.PG()
            log.debug('"%s" PG: %f' % (ingredient, pg))
            totalPg += pg
        return totalPg

    def gravity(self):
        return 1 + (self.PG() * PPG2HWE / (1000 * self.recipeVolume))

    def boilGravity(self):
        return 1 + (self.PG() * PPG2HWE / (1000 * self.boilVolume))

if __name__ == "__main__":
    recipe = Recipe()
##    g1 = Grain('grain', 37, 10)
##    i1 = Ingredient(g1, 1, 60)
##    recipe.addIngredient(i1)

    g2 = Grain("DME", potential = 45, colour = 12.5)
    i2 = Ingredient(g2, amount = 2.5)
    recipe.addIngredient(i2)

    h = Hops('bramling cross', 0.065)
    i = Ingredient(h, 0.05, 60)
    recipe.addIngredient(i)

    h = Hops('saaz', 0.035)
    i = Ingredient(h, 0.025, 30)
    recipe.addIngredient(i)
    #recipe.boilGravity = 1.041

    print 'gravity', recipe.gravity()
    print 'ibu', recipe.IBU()
