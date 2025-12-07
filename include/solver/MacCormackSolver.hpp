#ifndef MACCORMACKSOLVER_HPP
#define MACCORMACKSOLVER_HPP

#include <memory>

#include "Solver.hpp"
#include "bc/BoundaryManager.hpp"
#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"

/**
 * @class MacCormackSolver
 * @brief Two-step predictor–corrector (MacCormack) solver for the 1D Euler equations.
 *
 * This solver implements the classical explicit MacCormack scheme in conservative form
 * on a uniform 1D grid:
 *
 * Predictor step (forward difference):
 * \f[
 *   U_j^{*} = U_j^n
 *   - \frac{\Delta t}{\Delta x}\,\Big(F(U_{j+1}^n) - F(U_j^n)\Big),
 * \f]
 *
 * Corrector step (backward difference):
 * \f[
 *   U_j^{n+1} = \frac12\Bigg(
 *      U_j^n + U_j^{*}
 *      - \frac{\Delta t}{\Delta x}\,\Big(F(U_j^{*}) - F(U_{j-1}^{*})\Big)
 *   \Bigg),
 * \f]
 *
 * где \f$F(U)\f$ — физический поток Эйлера, вычисляемый через Primitive и EulerFlux.
 *
 * Используются:
 *  - DataLayer для хранения примитивных и вспомогательных полей,
 *  - EOS для преобразований Prim <-> Cons и давления,
 *  - TimeStepCalculator для CFL–шага по времени,
 *  - PositivityLimiter как «страховка» на каждом полном шаге.
 *
 * Схема:
 *  - консервативная по форме (разностный оператор строится из разностей потоков),
 *  - номинально второго порядка по времени и пространству для гладких решений,
 *  - достаточно простой и быстрой для тестов типа Sod.
 */
class MacCormackSolver : public Solver {
public:
    /**
     * @brief Конструктор по глобальным настройкам.
     *
     * Использует:
     *  - settings.gamma     — показатель адиабаты,
     *  - settings.cfl       — число CFL,
     *  - settings.t_end     — конечное время,
     *  - settings.N, L_x    — для вычисления шага dx (если заданы),
     *  - settings.dim       — размерность (ожидается dim = 1).
     *
     * @param settings Глобальные настройки задачи (сохраняются по значению).
     */
    explicit MacCormackSolver(const Settings& settings);

    /**
     * @brief Выполняет один шаг схемы МакКормака.
     *
     * Алгоритм:
     *  1. Применить граничные условия (BoundaryManager).
     *  2. Вычислить dx (по Settings либо по координатам xc в DataLayer).
     *  3. Посчитать шаг по времени dt из CFL через TimeStepCalculator.
     *  4. При необходимости «подрезать» dt, чтобы не превысить t_end.
     *  5. Выполнить предиктор–корректор для всех core–ячеек:
     *      - сформировать массивы консервативных переменных U^n,
     *      - посчитать физические потоки F(U^n),
     *      - получить предсказанные U^* (forward difference),
     *      - посчитать потоки F(U^*),
     *      - получить U^{n+1} по формуле корректировки (backward difference),
     *      - применить PositivityLimiter и сохранить результат в DataLayer.
     *  6. Увеличить t_cur на dt и вернуть dt.
     *
     * Если core-ячеек меньше двух или dt <= 0, шаг не делается и возвращается 0.0.
     *
     * @param layer DataLayer c текущим состоянием (модифицируется на месте).
     * @param t_cur Текущее время моделирования (увеличивается на dt).
     * @return Фактический шаг по времени dt, либо 0.0 если шаг не выполнен.
     */
    auto Step(DataLayer& layer, double& t_cur) -> double override;

    /**
     * @brief Устанавливает число CFL.
     *
     * @param cfl Новое значение CFL (ожидается 0 < cfl <= 1).
     */
    void SetCfl(double cfl) override { cfl_ = cfl; }

    /**
     * @brief Назначает граничные условия по оси.
     *
     * Для 1D задач ожидается axis = 0.
     *
     * @param axis     Индекс пространственной оси.
     * @param left_bc  ГУ на левой границе.
     * @param right_bc ГУ на правой границе.
     */
    void AddBoundary(int axis, std::shared_ptr<BoundaryCondition> left_bc,
                     std::shared_ptr<BoundaryCondition> right_bc) override;

private:
    /// Локальная копия настроек.
    Settings settings_;

    /// Менеджер граничных условий.
    BoundaryManager boundary_manager_;

    /// Минимально допустимая плотность (для PositivityLimiter).
    double rho_min_;

    /// Минимально допустимое давление (для PositivityLimiter).
    double p_min_;

    /**
     * @brief Вычисляет шаг по пространству dx для равномерной 1D сетки.
     *
     * Приоритет:
     *  1. Если settings_.N > 0 и settings_.L_x > 0: dx = L_x / N.
     *  2. Иначе попытаться взять dx из координат cell-centers:
     *       dx = xc[core_start+1] - xc[core_start].
     *  3. Если ничего не получилось, вернуть dx = 1.0.
     *
     * @param layer DataLayer, содержащий координаты xc.
     * @return Шаг по пространству dx.
     */
    [[nodiscard]] auto ComputeDx(const DataLayer& layer) const -> double;

    /**
     * @brief Записывает консервативное состояние в DataLayer.
     *
     * Переводит (rho, rhoU, E) в примитивные величины и вспомогательные поля:
     *  - rho, u, P
     *  - p  (импульс)
     *  - V  (удельный объем)
     *  - U  (удельная внутренняя энергия)
     *  - e  (полная энергия на объём)
     *  - m  (масса ячейки = rho * dx)
     *
     * Давление вычисляется через EOS::Pressure().
     *
     * @param uc    Консервативное состояние.
     * @param i     Индекс ячейки в линейной разметке (включая ghost).
     * @param dx    Шаг по пространству (для массы).
     * @param layer DataLayer, куда записать результат.
     */
    void StoreConservativeCell(const Conservative& uc, int i, double dx,
                               DataLayer& layer) const;
};

#endif  // MACCORMACKSOLVER_HPP
