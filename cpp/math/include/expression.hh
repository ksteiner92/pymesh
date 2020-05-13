//
// Created by klaus on 2019-12-06.
//

#ifndef PYULB_EXPRESSION_HH
#define PYULB_EXPRESSION_HH

#include <type_traits>

#include "utils.hh"

namespace Expr
{

template<typename T>
using is_complex = std::disjunction<std::is_same<std::remove_reference_t<T>, std::complex<float>>,
        std::is_same<std::remove_reference_t<T>, std::complex<double>>,
        std::is_same<std::remove_reference_t<T>, std::complex<long double>>>;

template<typename T>
static constexpr bool is_complex_v = is_complex<T>::value;

template<typename T>
using is_number = std::disjunction<std::is_arithmetic<std::remove_reference_t<T>>, is_complex<T>>;

template<typename T>
static constexpr bool is_number_v = is_number<T>::value;

template<typename T>
struct Expression;

template<typename Left, typename Right, typename Result = decltype(std::declval<Left>() * std::declval<Right>())>
struct Multiply;

template<typename Left, typename Right, typename Result = decltype(std::declval<Left>() / std::declval<Right>())>
struct Divide;

template<typename Left, typename Right, typename Result = decltype(std::declval<Left>() + std::declval<Right>())>
struct Add;

template<typename Left, typename Right, typename Result = decltype(std::declval<Left>() - std::declval<Right>())>
struct Subtract;

template<typename T>
struct Negate;

template<typename T>
struct Constant;

template<typename T>
struct Variable;

template<typename T, typename R = std::remove_cv_t<std::remove_reference_t<T>>>
std::enable_if_t<is_number_v<R>, Constant<R>>
expr(T&& number)
{
   return Constant<R>(std::move(number));
}

/*template<typename T, typename R = std::remove_reference_t<T>>
std::enable_if_t<std::is_base_of_v<Expression<R>, R>, R>
expr(T&& number)
{
   return number;
}*/

template<typename T>
auto eval(const Expression<T>& expr);

template<typename T>
struct Expression
{
   virtual inline T operator()(std::size_t i) const = 0;

   virtual std::size_t size() const = 0;

   /*template<typename Right, typename R = std::remove_reference_t<Right>>
   Multiply<T, R> operator*(Right&& r)
   {
      return Multiply<T, R>(std::move(*this), std::move(expr(r)));
   }*/

   template<typename Right>
   Variable<Right> operator=(const Expression<Right>& r)
   {
      return eval(r);
   }

   template<typename Right>
   Multiply<T, Right> operator*(const Expression<Right>& r)
   {
      return Multiply<T, Right>(std::move(*this), std::move(const_cast<Expression<Right>&>(r)));
   }

   template<typename Right>
   Divide<T, Right> operator/(const Expression<Right>& r)
   {
      return Divide<T, Right>(std::move(*this), std::move(const_cast<Expression<Right>&>(r)));
   }

   template<typename Right>
   Add<T, Right> operator+(const Expression<Right>& r)
   {
      return Add<T, Right>(std::move(*this), std::move(const_cast<Expression<Right>&>(r)));
   }

   template<typename Right>
   Subtract<T, Right> operator-(const Expression<Right>& r)
   {
      return Subtract<T, Right>(std::move(*this), std::move(const_cast<Expression<Right>&>(r)));
   }

   template<typename Right>
   Negate<T> operator-()
   {
      return Negate<T>(std::move(*this));
   }
};

template<typename Left, typename Right, typename Result>
struct BinaryOperator : public Expression<Result>
{
   BinaryOperator(Expression<Left>&& l, Expression<Right>&& r) : l(std::forward<Expression<Left>&>(l)), r(std::forward<Expression<Right>&>(r)) {}

   std::size_t size() const override
   {
      return std::max(l.size(), r.size());
   }

   Expression<Left>& l;
   Expression<Right>& r;
};

template<typename T>
struct UnaryOperator : public Expression<T>
{
   explicit UnaryOperator(Expression<T>&& expr) : expr(expr){}

   std::size_t size() const override
   {
      return expr.size();
   }

   Expression<T>& expr;
};

template<typename T>
struct Negate : public UnaryOperator<T>
{
   using super = UnaryOperator<T>;
   using super::UnaryOperator;

   T operator()(std::size_t i) const override
   {
      return -super::expr(i);
   }
};

template<typename Left, typename Right, typename Result>
struct Add : public BinaryOperator<Left, Right, Result>
{
   using super = BinaryOperator<Left, Right, Result>;
   using super::BinaryOperator;

   Result operator()(std::size_t i) const override
   {
      return super::l(i) + super::r(i);
   }
};

template<typename Left, typename Right, typename Result>
struct Subtract : public BinaryOperator<Left, Right, Result>
{
   using super = BinaryOperator<Left, Right, Result>;
   using super::BinaryOperator;

   Result operator()(std::size_t i) const override
   {
      return super::l(i) - super::r(i);
   }
};

template<typename Left, typename Right, typename Result>
struct Multiply : public BinaryOperator<Left, Right, Result>
{
   using super = BinaryOperator<Left, Right, Result>;
   using super::BinaryOperator;

   Result operator()(std::size_t i) const override
   {
      return super::l(i) * super::r(i);
   }
};

template<typename Left, typename Right, typename Result>
struct Divide : public BinaryOperator<Left, Right, Result>
{
   using super = BinaryOperator<Left, Right, Result>;
   using super::BinaryOperator;

   Result operator()(std::size_t i) const override
   {
      return super::l(i) / super::r(i);
   }
};

template<typename Expr, template<typename> class Result>
using enable_if_expression =
std::enable_if_t<std::is_base_of_v<Expression, std::remove_reference_t<Expr>>,
                 Result<std::remove_reference_t<Expr>>>;

template<typename Left, template<typename, typename> class Operator, typename Right, typename Result = decltype(std::declval<Operator>()(0))>
using enable_if_binary_expression =
std::enable_if_t<std::conjunction_v<std::is_base_of<Expression<Result>, std::remove_reference_t<Left>>,
        std::is_base_of<Expression<Result>, std::remove_reference_t<Right>>,
        std::is_base_of<BinaryOperator<std::remove_reference_t<Left>, std::remove_reference_t<Right>, Result>,
                              Operator<std::remove_reference_t<Left>, std::remove_reference_t<Right>>>>,
        Operator<std::remove_reference_t<Left>, std::remove_reference_t<Right>>>;

template<typename T>
struct Constant : public Expression<T>
{
   explicit Constant(T&& value) : value(std::forward<T>(value)) {}

   T operator()(std::size_t i) const override
   {
      return value;
   }

   std::size_t size() const override
   {
      return 0;
   }

   T value;
};

template<typename T>
struct Variable : public Expression<T>
{
   explicit Variable(std::valarray<T>&& value) : value(std::forward<std::valarray<T>>(value)) {}

   explicit Variable(std::size_t s) : value(std::valarray<T>(s)) {}

   T operator()(std::size_t i) const override
   {
      return value[i];
   }

   const std::valarray<T>& operator()() const
   {
      return value;
   }

   std::size_t size() const override
   {
      return value.size();
   }


   T& operator[](std::size_t i)
   {
      return value[i];
   }

   std::valarray<T> value;
};

/*template<typename Left, template<typename, typename> class Operator, typename Right, typename Result = decltype(std::declval<Operator>()(0))>
struct BinaryExpression : public Expression<Result>
{
   BinaryExpression(Left&& l, Right&& r) : op(std::forward<Left>(l), std::forward<Right>(r)) {}

   Result operator()(std::size_t i) const override
   {
      return op(i);
   }

   Operator<Left, Right> op;
};*/

template<typename T>
auto eval(const Expression<T>& expr)
{
   Variable<T> res(expr.size());
   for (std::size_t i = 0; i < expr.size(); ++i) {
      res[i] = expr(i);
   }
   return res;
}

/*template<typename Left, typename Right, typename L = std::remove_reference_t<Left>, typename R = std::remove_reference_t<Right>>
enable_if_binary_expression<L, Multiply, R>
operator*(Left&& l, Right&& r)
{
   return enable_if_binary_expression<L, Multiply, R>(std::move(l), std::move(r));
}

template<typename Left, typename Right, typename L = std::remove_reference_t<Left>, typename R = std::remove_reference_t<Right>>
enable_if_binary_expression<L, typename Add<L, R>, R>
operator+(Left&& l, Right&& r)
{
   return enable_if_binary_expression<L, Add, R>(std::move(l), std::move(r));
}

template<typename Expr>
enable_if_expression<Expr, Negate>
operator-(Expr&& expr)
{
   return enable_if_expression<Expr, Negate>(std::move(expr));
}

template<typename Left, typename Right>
auto operator-(Left&& l, Right&& r)
{
   return l + (-r);
}*/

/*template<typename Expr>
enable_if_binary_expression<Expr, Multiply, Expr>
sq(Expr&& expr)
{
   return enable_if_binary_expression<Expr, Multiply, Expr>(std::move(expr), std::move(expr));
}*/

}

#endif //PYULB_EXPRESSION_HH
