#ifndef est_for_nde_hpp
#define est_for_nde_hpp

#include <cstdio>
#include <iterator>
#include <vector>

static const int intNULL = 2111111111;

class IntIterator : public std::iterator<std::bidirectional_iterator_tag, int, void>
{
public:
  IntIterator(int* begin)
    : m_pos(begin)
  {};

  IntIterator(int* begin, int)
    : m_pos(begin)
  {};

  IntIterator()
    : m_pos(NULL)
  {};

  int& operator *() const {
    return *m_pos;
  };

  IntIterator& operator ++(){
    ++m_pos;
    return *this;
  };

  IntIterator operator ++(int){
    int* p = m_pos;
    ++(*this);
    return IntIterator(p);
  };

  IntIterator& operator --(){
    --m_pos;
    return *this;
  };

  const bool operator !=(const IntIterator& ite) const {
    if(ite.m_pos == NULL && *m_pos == intNULL) return false;
    return (m_pos != ite.m_pos);
  };

  const bool operator ==(const IntIterator& ite) const {
    if(ite.m_pos == NULL && *m_pos == intNULL) return true;
    return (m_pos == ite.m_pos);
  };

private:
  int* m_pos;
};
  

#define for_nde(i,abcABC) \
if(int jabcABC=0 == 0)for(IntIterator i((int *)&abcABC);jabcABC++<7; ++i)

#define for_nde2(i,abcABC) \
if(int jabcABC2=0 == 0)for(IntIterator i((int *)&abcABC);jabcABC2++<7; ++i)

#define for_abc(i,abcABC) \
if(int jabc=0 == 0)for(IntIterator i((int *)&abcABC);jabc++<4; ++i)

#define for_ABC(i,abcABC) \
if(int jABC=0 == 0)for(IntIterator i((int *)&abcABC);jABC++<7; ++i)if(4<jABC)

#define for_ABC2(i,abcABC) \
if(int jABC2=0 == 0)for(IntIterator i((int *)&abcABC);jABC2++<7; ++i)if(4<jABC2)


#endif
