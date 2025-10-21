#ifndef _MES_CORE_H
#define _MES_CORE_H
#define INTEGRATION_POINTS 2


namespace MES {
	typedef double real_t;

	struct node {
		union {
			struct {real_t x, y;};
			struct {real_t ksi, eta;};
		};
	};
}
#endif // _MES_CORE_H