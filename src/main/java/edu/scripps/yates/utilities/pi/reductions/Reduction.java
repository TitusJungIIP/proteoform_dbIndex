/*
 *  Copyright (C) 2009 Nasser Giacaman, Oliver Sinnen
 *
 *  This file is part of Parallel Iterator.
 *
 *  Parallel Iterator is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or (at 
 *  your option) any later version.
 *
 *  Parallel Iterator is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
 *  Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along 
 *  with Parallel Iterator. If not, see <http://www.gnu.org/licenses/>.
 */

package edu.scripps.yates.utilities.pi.reductions;

/**
 * Defines a reduction and includes a range of built-in reductions (only a few
 * common ones are implemented). Users may also define their own by implementing
 * this interface.
 *
 * @author Nasser Giacaman
 * @author Oliver Sinnen
 */
public interface Reduction<E> {

	/**
	 * Specifies a reduction as defined by 2 elements into 1.
	 * 
	 * The reduction must obey the following 2 constraints:
	 * <ul>
	 * <li><i>Associative</i>: the order of evaluating the reduction makes no
	 * difference, and
	 * <li><i>Commutative</i>: the order of the thread-local values makes no
	 * difference.
	 * </ul>
	 * 
	 * @param first
	 *            The first element in the reduction.
	 * @param second
	 *            The second element in the reduction.
	 * @return The result of reducing <code>first</code> with
	 *         <code>second</code>.
	 */
	public E reduce(E first, E second);

	/*
	 * TODO A few reductions are implemented below. More should be added
	 */

	/**
	 * Returns the maximum for <code>Integer</code>
	 */
	public static Reduction<Integer> IntegerMAX = new Reduction<Integer>() {
		@Override
		public Integer reduce(Integer first, Integer second) {
			return Math.max(first, second);
		}
	};

	/**
	 * Returns the minimum for <code>Integer</code>
	 */
	public static Reduction<Integer> IntegerMIN = new Reduction<Integer>() {
		@Override
		public Integer reduce(Integer first, Integer second) {
			return Math.min(first, second);
		}
	};

	/**
	 * Returns the sum for <code>Integer</code>
	 */
	public static Reduction<Integer> IntegerSUM = new Reduction<Integer>() {
		@Override
		public Integer reduce(Integer first, Integer second) {
			return first + second;
		}
	};

	/**
	 * Returns the maximum for <code>Double</code>
	 */
	public static Reduction<Double> DoubleMAX = new Reduction<Double>() {
		@Override
		public Double reduce(Double first, Double second) {
			return Math.max(first, second);
		}
	};

	/**
	 * Returns the minimum for <code>Double</code>
	 */
	public static Reduction<Double> DoubleMIN = new Reduction<Double>() {
		@Override
		public Double reduce(Double first, Double second) {
			return Math.min(first, second);
		}
	};

	/**
	 * Returns the sum for <code>Double</code>
	 */
	public static Reduction<Double> DoubleSUM = new Reduction<Double>() {
		@Override
		public Double reduce(Double first, Double second) {
			return first + second;
		}
	};

}
