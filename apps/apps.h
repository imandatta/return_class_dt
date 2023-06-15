#pragma once

#include <iostream>
#include <memory>
#include <array>
#include <string>


#include "definitions.h"

namespace apps
{

class WriteoutDiagnostic
{
private:
    real _dt;
    std::array<real, 3> _x;
    std::string _physics;

public:
    // Constructor
    WriteoutDiagnostic()
    {
    }
    // Constructor
    WriteoutDiagnostic(real dt, const std::array<real, 3>& x, const std::string& physics)
        : _dt(dt), _x(x), _physics(physics)
    {
    }

    // Minimum dt
    static const WriteoutDiagnostic& minDt(const WriteoutDiagnostic& a,
                                           const WriteoutDiagnostic& b)
    {
        return a._dt < b._dt ? a : b;
    }

    // Setters
    void setDt(real dt)
    {
        _dt = dt;
    }
    void setX(const std::array<real, 3>& x)
    {
        _x = x;
    }
    void setPhysics(const std::string& physics)
    {
        _physics = physics;
    }

    void setValues(real dt, const std::array<real, 3>& x, const std::string& physics)
    {
        setDt(dt);
        setX(x);
        setPhysics(physics);
    }

    // Getters
    real getDt() const
    {
        return _dt;
    }

    const std::array<real, 3>& getX() const
    {
        return _x;
    }

    const std::string& getPhysics() const
    {
        return _physics;
    }
};

class AppBase
{
protected:
    std::shared_ptr<WriteoutDiagnostic> _wd;
    std::string _physics;

public:
    // Constructor
    AppBase()
    {
        _wd = std::make_shared<WriteoutDiagnostic>();
    }

    virtual void setup()
    {
    }

    virtual const std::shared_ptr<WriteoutDiagnostic>
    call(const std::array<real, 3>& x) const
    {
        std::cerr << "Unimplemented call()\n";
        exit(-1);
    }
};

class App0 : public AppBase
{
private:
public:
    // Constructor
    App0()
    {
        _physics = "app0";
    }

    // Setup
    void setup() override
    {
        AppBase::setup();
    }

    const std::shared_ptr<WriteoutDiagnostic>
    call(const std::array<real, 3>& x) const override;
};

class App1 : public AppBase
{
private:
public:
    // Constructor
    App1()
    {
        _physics = "app1";
    }

    // Setup
    void setup() override
    {
        AppBase::setup();
    }

    const std::shared_ptr<WriteoutDiagnostic>
    call(const std::array<real, 3>& x) const override;
};

class App2 : public AppBase
{
private:
public:
    // Constructor
    App2()
    {
        _physics = "app2";
    }

    // Setup
    void setup() override
    {
        AppBase::setup();
    }

    const std::shared_ptr<WriteoutDiagnostic>
    call(const std::array<real, 3>& x) const override;
};

} // namespace apps
